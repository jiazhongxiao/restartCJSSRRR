function [eigenvalues, eigenvectors, flag, MVs, info] = CJSSRRR_hybrid(A, interval, varargin)
%   CJSSRRR_hybrid: A refined CJ--SS--RR algorithm with selectable restarts
%
%   input param:
%       A:          n-by-n Hermitian matrix, the matrix to be solved
%       interval:   [a, b], the interval of interest
%       varargin:   Optional parameters, passed as name-value pairs or a structure
%   name-value pair:
%       moment:     Number of moments                     [integer   :     {8}]
%       maxItr:     maximum number of iterations          [integer   :    {50}]
%       eigTol:     Residual tolerance in [0,1)           [double    : {1e-12}]
%       eigenspace: Desired eigenspace (for testing)      [matrix    :    {[]}]
%       spectrum:   Spectrum of A                         [vector    :    {[]}]
%       nev:        Number of eigenvalues in the region   [integer   :     {0}]
%       V:          Starting matrix to generate S         [matrix    :    {[]}]
%       degree:     Degree of the CJ polynomial series    [integer   :     {0}]
%       degreeD:    Positive D (recommended 1~10)         [double    :     {1}]
%       omega:      Scaling factor for the dimension of
%                     search subspace, dim = nev * omega  [double    :   {1.5}]
%       disp:       Diagnostic information display level  [integer   :     {0}]
%       restart:    Restart format                        ["natural" | "combine" | {"hybrid"}]
%       orthMethod: Search-space orthonormalization       [  {"QR"}  |     "TSVD" ]
%       svdTol:     Relative TSVD truncation threshold    [double    :   {100*eps}]
%       eps_switch: Hybrid switching threshold in [0,1)   [double    : {eigTol^.5}]
%
%   output param:
%       eigenvalues:             Retained refined Ritz values within the interval
%       eigenvectors:            Corresponding refined Ritz vectors
%       flag:                    Algorithm exit flag, -1: natural restart residual growth,
%                                  0: converged, 1: maxItr reached, 2: no iterations
%       MVs:                     Total matrix-vector products count
%       info:                    Structure containing additional information
%       ____.time:               Timing information for various stages
%       _________.spectrum:      Time for spectrum estimation
%       _________.CJcoeffs:      Time for CJ coefficients computation
%       _________.nev:           Time for nev estimation
%       _________.subspace:      Time for subspace generation, [MV-time flops-time]
%       _________.orth:          Time for full search-space orthonormalization
%       _________.RR:            Time for Rayleigh-Ritz
%       _________.refined:       Time for refined procedure, [preprocess-time compute-time]
%       _________.residual:      Time for residual computation
%       _________.total:         Total time for the algorithm
%       ____.spectrum:           Estimated spectrum of A
%       ____.degree:             Degree of the CJ polynomial series
%       ____.degreeD:            D used in degree estimation
%       ____.nev:                Estimated number of eigenvalues in the region
%       ____.V:                  Starting matrix V used for subspace generation
%       ____.moment:             Number of moments
%       ____.dimension:          Number of columns in the initial raw search matrix S
%       ____.iteration:          Number of iterations performed
%       ____.restart:            Normalized requested restart format
%       ____.orthMethod:         Effective search-space orthonormalization method
%       ____.svdTol:             Relative TSVD truncation threshold
%       ____.eps_switch:         Hybrid switching threshold
%       ____.switch:             Hybrid switch information, [itr, cnd]
%       ____.itrInfo:            Structure containing iteration information
%       ____________.dimRR:      Rayleigh-Ritz projection dimension
%       ____________.residual:   Maximum relative residual norms
%       ____________.nev:        Number of retained refined Ritz values
%       ____________.condition:  Condition numbers of retained vectors
%       ____________.insideRes:  Residuals of approximations inside the interval
%       ____________.insideNum:  Number of approximations inside the interval
%       ____________.sinSearch:  Sine of angles for the pre-truncation search subspace
%                                  at each iteration
%       ____________.sinRetain:  Sine of angles for retained refined Ritz vectors

t0 = tic; t1  = t0;
MVs_t = 0; % The count of matrix-vector products
flops = 0; % The count of other flops of CJ approximation

%%  Check input parameters
[moment, maxItr, eigTol, eigenspace, spectrum, nev, V, ...
    degree, degreeD, omega, disp, restart, orthMethod, ...
    svdTol, eps_switch] = CheckInputs(A, interval, varargin{:});
isGivenSpectrum   = ~isempty(spectrum);
isGivenNev        = nev > 0;
isGivenV          = ~isempty(V);
isGivenEigenspace = ~isempty(eigenspace);


%%  Default parameters
maxSpectrumTry  = 5;  % Maximum number of attempts to estimate the spectrum
SpectrumEstStep = 20; % Default step for spectrum estimation
NevEstimateSize = 20; % Default size of the random matrix for nev estimation


%%  Prepare internal and output variables
colA = size(A, 2); % Number of columns in A
info = struct();   % Initialize info structure
info.time = struct("spectrum",     0, "CJcoeffs", 0,   "nev", 0, ...
                   "subspace", [0 0],     "orth", 0,    "RR", 0, ...
                   "refined",  [0 0], "residual", 0, "total", 0);
info.restart    = restart;
info.orthMethod = orthMethod;
info.svdTol     = svdTol;
info.eps_switch = eps_switch;
info.switch     = struct("itr", 0, "cnd", 0);
info.iteration  = 0;


%   Estimate the spectrum of A until overestimation
for estItr = 1:maxSpectrumTry
%%  Estimate the spectrum of A and map it to [-1, 1]
    time_spectrum = tic;
    % Estimate the spectrum of A
    if ~isGivenSpectrum
        spectrum = EstimateSpectrum(A, SpectrumEstStep);
        MVs_t = MVs_t + SpectrumEstStep;
    end
    info.time.spectrum = info.time.spectrum + toc(time_spectrum);
    % Compute the center and radius of the spectrum
    center = (min(spectrum) + max(spectrum)) / 2;
    radius = (max(spectrum) - min(spectrum)) / 2;
    % Map the spectrum to [-1, 1]
    % A_scaled = (A - center * speye(size(A))) / radius;
    % Update the interval of interest
    a = (interval(1) - center) / radius;
    b = (interval(2) - center) / radius;
    % Ensure the updated interval is within [-1, 1]
    a = max(-1, min(1, a));
    b = max(-1, min(1, b));

%%  Estimate the degree of the CJ series and compute the CJ coefficients
    % Estimate the degree of the CJ polynomial series
    degree = EstimateDegree(a, b, moment, degree, degreeD);
    % Compute the CJ coefficients
    time_CJcoeffs = tic;
    coeffs = GetCJCoeffs(a, b, degree, moment);
    info.time.CJcoeffs = info.time.CJcoeffs + toc(time_CJcoeffs);

%%  Test whether the spectrum is well estimated
    if isGivenSpectrum; break; end
    v = randn(colA, 1);
    [nev_est, ~, ~, info.time, MVs_t, flops] = ...
        EstimateNev(A, center, radius, v, coeffs(:,1), info.time, MVs_t, flops);
    clear v;
    if abs(nev_est) < colA; break; end
    if estItr == maxSpectrumTry
        error("CJSSRRR_hybrid:EstimateSpectrumFailed", ...
            "Failed to estimate the spectrum of A after %d attempts. ", ...
            maxSpectrumTry);
    end
end
%   Record the spectrum and degree in info
info.spectrum = spectrum;
info.degree   = degree;
info.degreeD  = degreeD;


%%  Estimate the number of eigenvalues in the region
if ~isGivenNev % If nev is not given, estimate it
    if ~isGivenV
        V = randn(colA, NevEstimateSize); % Generate a random matrix V
    end
    time_nev = tic;
    [nev, S, V, info.time, MVs_t, flops] = ...
        EstimateNev(A, center, radius, V, coeffs, info.time, MVs_t, flops);
    info.time.nev = info.time.nev + toc(time_nev);
end
%   Record the estimated nev in info
info.nev = nev;


%%  Exit if maxItr is 0
if maxItr <= 0
    eigenvalues  = [];
    eigenvectors = [];
    flag         = 2; % No iterations were performed
    MVs = [MVs_t flops/(2*nnz(A))];
    info.time.total = toc(t0);
    return;
end


%%  Resize the dimension of search subspace for the first iteration
if ~isGivenNev
    if ~isGivenV
        colV = ceil(nev * omega / moment); % Number of columns in V
        [S, V, info.time, MVs_t, flops] = ...
            ResizeCJSubspace(A, center, radius, S, V, colV, coeffs, info.time, MVs_t, flops);
    end
else % nev > 0
    if ~isGivenV
        colV = ceil(nev * omega / moment); % Number of columns in V
        V = randn(colA, colV);
    end
    [S, info.time.subspace, MVs_t, flops] = ...
        GetCJSubspace(A, center, radius, V, coeffs, info.time.subspace, MVs_t, flops);
end
%   Record the starting matrix V, the moment size and the dimension of S in info
info.V         = V;
info.moment    = moment;
info.dimension = size(S, 2);


%%  Prepare parameters for the main loop
a      = interval(1);                  % Left endpoint of the interval
b      = interval(2);                  % Right endpoint of the interval
colV   = size(V, 2);                   % Number of columns in V
nrmA   = max(abs(spectrum));           % Norm of A
minGap = max(eigTol, eps*1024) * nrmA; % Minimum gap for eigenvalues
runSI  = strcmpi(restart, "hybrid");   % Whether the hybrid method remains natural


%%  Initialize iteration information variables
vecDim = zeros(maxItr, 1); % Actual projection dimension in each iteration
vecRes =  ones(maxItr, 1); % Maximum relative residual norms
vecNev = zeros(maxItr, 1); % Estimated nev in each iteration
vecCnd = zeros(maxItr, 1); % Condition numbers of the desired Ritz vectors
res_in =  ones(maxItr, 1); % Maximum relative residual of inside approximations
num_in = zeros(maxItr, 1); % Number of inside approximations
% Sine of the angles between desired eigenspace and R(S) or R(X)
sine_S =  ones(maxItr, 1);
sine_X =  ones(maxItr, 1);


%%  Display initial information
if disp > 0
    fprintf("%1$s CJSSRRR_hybrid %1$s\n", sprintf("%s", repmat("-", 1, 20)));

    fprintf("CJSSRRR_hybrid: Spectrum: [%.2e, %.2e], nev: %.1f.\n", spectrum, nev);
    fprintf("\t Degree: %d, D: %.1f. Moment: %d.\n", ...
            degree, degreeD, moment);
    fprintf("\t Block column size: %d, search dim.: %d.\n", colV, size(S,2));
    fprintf("\t Restart: %s restart", restart);
    if strcmpi(restart, "hybrid")
        fprintf(", eps_switch=%.2e", eps_switch);
    end
    fprintf(".\n");
    fprintf("\t Orth.: %s", orthMethod);
    if strcmpi(orthMethod, "TSVD")
        if svdTol == 0
            fprintf(" (svdTol=0, no truncation)");
        else
            fprintf(" (svdTol=%.2e)", svdTol);
        end
    end
    fprintf(".\n");
end
if disp > 1
    fprintf(" #it | residual | nev | cond(X) | dim |  res_in  | n_in |");
    if disp > 2 && isGivenEigenspace; fprintf(" sin_S | sin_X |"); end
    fprintf(" time(min)\n");
end


%%%%%%%%%% Main loop for the CJSSRRR_hybrid iterations %%%%%%%%%%
for itr = 1:maxItr
    %% Compute the orthonormal basis of the search subspace
    if itr > 1
        [S, info.time.subspace, MVs_t, flops] = ...
            GetCJSubspace(A, center, radius, V, coeffs, info.time.subspace, MVs_t, flops);
    end
    if isempty(S) || ~any(S(:))
        error("CJSSRRR_hybrid:EmptySearchSubspace", ...
            "The search matrix S is empty or contains only zeros.");
    end
    time_orth = tic;
    switch orthMethod
        case "TSVD"
            [Q, sigma, ~] = svd(S, "econ", "vector");
            if svdTol > 0 && ~isempty(sigma)
                Q = Q(:, sigma >= norm(S,"fro")*svdTol);
            end
        case "QR"
            [Q, ~] = qr(S, "econ");
    end
    info.time.orth = info.time.orth + toc(time_orth);
    if isempty(Q) || ~any(Q(:))
        error("CJSSRRR_hybrid:EmptySearchSubspace", ...
            "The orthogonalized search matrix Q is empty or contains only zeros.");
    end
    vecDim(itr) = size(Q, 2);

    %% Compute the Ritz pairs and the residual norms
    % Rayleigh-Ritz procedure
    time_RR = tic;
    [Y, E, AQ, QAQ] = RayleighRitz(A, Q);
    MVs_t = MVs_t + size(Q, 2); % Count the matrix-vector products for A*Q
    info.time.RR = info.time.RR + toc(time_RR);
    % Residual norms
    time_res = tic;
    X   = Q * Y;                % Ritz vectors
    res = vecnorm(AQ*Y - X.*E); % Residual norms
    info.time.residual = info.time.residual + toc(time_res);

    %% Select the approximations in the interval
    eig_in = E >= a & E <= b;

    %% Record the information of search subspace
    if isGivenEigenspace
        % Record the angle between the desired eigenspace and the initial search subspace
        sine_S(itr) = SinSpace(S, eigenspace);
    end

    %% Refined procedure
    [Y, E, time_refined, res, idx_copy] = ...
        RefinedRitz(Q, AQ, QAQ, E(eig_in), res(eig_in), minGap, [a b]);
    info.time.refined = info.time.refined + time_refined;
    % Update the Ritz vectors and residual norms
    time_res = tic;
    X = Q * Y; % Refined Ritz vectors
    info.time.residual = info.time.residual + toc(time_res);
    % Reselect the approximations in the interval
    eig_in = E >= a & E <= b;
    % Select non-copy vectors
    idx_want  = ~idx_copy; % non-copy vectors

    %% Record the iteration information
    % Select desired approximations
    eig_desired = eig_in & idx_want;
    % Record the information of desired approximations
    if nnz(eig_desired) > 0
        vecRes(itr) = max(res(eig_desired)) / nrmA;
        vecNev(itr) = nnz(eig_desired);
        vecCnd(itr) = cond(X(:, eig_desired));
    end
    % Record the information about inside approximations
    if nnz(eig_in) > 0
        res_in(itr) = max(res(eig_in)) / nrmA;
        num_in(itr) = nnz(eig_in);
    end
    % Record the retained-subspace angle only when retained vectors exist
    if isGivenEigenspace && nnz(eig_desired) > 0
        sine_X(itr) = SinSpace(X(:,eig_desired), eigenspace);
    end

    %% Test the convergence
    if vecRes(itr) < eigTol; flag = 0; break; end % Converged
    if itr == maxItr;        flag = 1; break; end % Reached maximum iterations

    %% Natural restart residual-growth exit
    if strcmpi(restart, "natural") && eigTol > eps && ...
            vecRes(itr) > 10*min(vecRes(1:itr))
        flag = -1;
        break;
    end

    %% Hybrid restart state transition
    if strcmpi(restart, "hybrid") && runSI
        if vecRes(itr) <= eps_switch
            runSI = false;
            info.switch.itr = itr;
            info.switch.cnd = 1;
        elseif vecRes(itr) > min(vecRes(1:itr))
            runSI = false;
            info.switch.itr = itr;
            info.switch.cnd = 2;
        end
    end

    %% Restart
    if moment == 1
        V = Q;
        colV = size(V, 2);
    else
        if strcmpi(restart, "natural") || runSI
            [V, ~] = qr(S(:,1:colV), "econ");
        else
            V = X(:,eig_in) * randn(nnz(eig_in), colV);
        end
    end

    %% Display the iteration information
    if disp > 1
        fprintf(" %3d | %.2e | %3d | %.1e | %3d | %.2e |  %3d |", ...
                itr, vecRes(itr), vecNev(itr), vecCnd(itr), ...
                vecDim(itr), res_in(itr), num_in(itr));
        if disp > 2 && isGivenEigenspace
            fprintf(" %.0e | %.0e |", sine_S(itr), sine_X(itr));
        end
        fprintf("  %.2f\n", toc(t1)/60); t1 = tic;
    end
end


%%  Finalize the output parameters
%   Desired eigenpairs
eigenvalues  = E(  eig_desired);
eigenvectors = X(:,eig_desired);
%   Total matrix-vector products count
MVs = [MVs_t flops/(2*nnz(A))];
%   Record the number of iterations
info.iteration = itr;
%   Record the information about the iterations
info.itrInfo = struct("residual",  vecRes(1:itr), "nev",       vecNev(1:itr), ...
                      "condition", vecCnd(1:itr), "insideRes", res_in(1:itr), ...
                      "insideNum", num_in(1:itr), "sinSearch", sine_S(1:itr), ...
                      "sinRetain", sine_X(1:itr), "dimRR",     vecDim(1:itr));
%   Record the time information
info.time.total = toc(t0);

%%  Display the result
if disp > 1
    fprintf(" %3d | %.2e | %3d | %.1e | %3d | %.2e |  %3d |", ...
            itr, vecRes(itr), vecNev(itr), vecCnd(itr), ...
            vecDim(itr), res_in(itr), num_in(itr));
    if disp > 2 && isGivenEigenspace
        fprintf(" %.0e | %.0e |", sine_S(itr), sine_X(itr));
    end
    fprintf("  %.2f\n", toc(t1)/60);
end
if disp > 0
    digit = ceil(log10(sum(MVs)));
    fprintf("Result: \n\tnev: %d, iteration: %d\n", vecNev(itr), itr);
    fprintf("\tResidual: %.7e\n", vecRes(itr));
    fprintf(sprintf("\tMVs cost:  %%%dd\n",digit), MVs(1));
    fprintf(sprintf("\tFlops/MV:  %%%d.2f\n",digit+3), MVs(2));
    fprintf("\tTotal MVs: %.2f\n", sum(MVs));
    fprintf("\tTotal time: %.4f min\n", info.time.total/60);
end

end %~ end of CJSSRRR_hybrid



%%  Subroutines
%   CheckInputs: Check and parse input parameters
function [moment, maxItr, eigTol, eigenspace, spectrum, nev, V, ...
        degree, degreeD, omega, disp, restart, orthMethod, ...
        svdTol, eps_switch] = CheckInputs(A, interval, varargin)
%   Indispensable parameters
%       moment:     (can't est.) Number of moments                     [integer   :     {8}]
%       maxItr:     (can't est.) Maximum number of iterations          [integer   :    {50}]
%       eigTol:     (can't est.) Residual tolerance in [0,1)           [double    : {1e-12}]
%   Optional parameters
%       eigenspace: (can't est.) Desired eigenspace (for testing)      [matrix    :    {[]}]
%       spectrum:                Spectrum of A                         [vector    :    {[]}]
%       nev:                     Number of eigenvalues in the region   [integer   :     {0}]
%       V:                       Starting matrix to generate S         [matrix    :    {[]}]
%       degree:                  Degree of the CJ polynomial series    [integer   :     {0}]
%       degreeD:                 Positive D (recommended 1~10)         [double    :     {1}]
%       omega:      (can't est.) Scaling factor for the dimension of
%                                search subspace, dim = nev * omega    [double    :   {1.5}]
%       disp:       (can't form) Diagnostic information display level  [integer   :     {0}]
%       restart:    (can't form) Restart format                        ["natural" | "combine" | {"hybrid"}]
%       orthMethod: (can't form) Search-space orthonormalization       [  {"QR"}  |     "TSVD" ]
%       svdTol:                  Relative TSVD truncation threshold    [double    :   {100*eps}]
%       eps_switch:              Hybrid switching threshold in [0,1)   [double    : {eigTol^.5}]

    % Check the size of the input matrix A
    if ~ishermitian(A)
        error("CJSSRRR_hybrid:InvalidMatrix", ...
            "Input matrix A must be Hermitian.");
    end
    colA = size(A, 2);

    % Check the interval
    if ~isvector(interval) || length(interval) ~= 2 || interval(1) >= interval(2)
        error("CJSSRRR_hybrid:InvalidInterval", ...
            "Interval must be a vector of length 2 with the first element less than the second.");
    end

    % Build input parser
    parser = inputParser;
    parser.KeepUnmatched   = true;  % Keep unmatched parameters
    parser.PartialMatching = false; % Do not accept unmatched abbreviations
    parser.StructExpand    = true;  % Use native structure input parsing

    % Default values
    defaultMoment     = 8;
    defaultMaxItr     = 50;
    defaultEigTol     = 1e-12;
    defaultEigenspace = [];
    defaultSpectrum   = [];
    defaultNev        = 0;
    defaultV          = [];
    defaultDegree     = 0;
    defaultDegreeD    = 1;
    defaultOmega      = 1.5;
    defaultDisp       = 0;
    defaultRestart    = "hybrid";
    defaultOrthMethod = "QR";
    defaultSvdTol     = 100*eps;
    % Parser placeholder; the effective default is sqrt(eigTol)
    defaultEps_switch = 0;

    % Define test functions
    isNumber    = @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x);
    isInteger   = @(x) isNumber(x) && x == fix(x);
    isCharStr   = @(x) (ischar(x) && isrow(x)) || (isstring(x) && isscalar(x));
    isMatrix    = @(x) isnumeric(x) && ismatrix(x);
    isVector    = @(x) isnumeric(x) && (isempty(x) || (isvector(x) && isreal(x)));
    isUnitRange = @(x) x >= 0 && x < 1;

    % Add parameters to the parser
    parser.addParameter("moment",     defaultMoment,      isInteger);
    parser.addParameter("maxItr",     defaultMaxItr,      isInteger);
    parser.addParameter("eigTol",     defaultEigTol,       isNumber);
    parser.addParameter("eigenspace", defaultEigenspace,   isMatrix);
    parser.addParameter("spectrum",   defaultSpectrum,     isVector);
    parser.addParameter("nev",        defaultNev,         isInteger);
    parser.addParameter("V",          defaultV,            isMatrix);
    parser.addParameter("degree",     defaultDegree,      isInteger);
    parser.addParameter("degreeD",    defaultDegreeD,      isNumber);
    parser.addParameter("omega",      defaultOmega,        isNumber);
    parser.addParameter("disp",       defaultDisp,        isInteger);
    parser.addParameter("restart",    defaultRestart,     isCharStr);
    parser.addParameter("orthMethod", defaultOrthMethod,  isCharStr);
    parser.addParameter("svdTol",     defaultSvdTol,       isNumber);
    parser.addParameter("eps_switch", defaultEps_switch,   isNumber);

    % Parse the input parameters
    parser.parse(varargin{:});

    % Extract the parameters
    moment     = parser.Results.moment;
    maxItr     = parser.Results.maxItr;
    eigTol     = parser.Results.eigTol;
    eigenspace = parser.Results.eigenspace;
    spectrum   = parser.Results.spectrum;
    nev        = parser.Results.nev;
    V          = parser.Results.V;
    degree     = parser.Results.degree;
    degreeD    = parser.Results.degreeD;
    omega      = parser.Results.omega;
    disp       = parser.Results.disp;
    restart    = parser.Results.restart;
    orthMethod = parser.Results.orthMethod;
    svdTol     = parser.Results.svdTol;
    eps_switch = parser.Results.eps_switch;

    % Validate the parameters
    % Check moment
    moment = UseDefault(moment, moment > 0, defaultMoment, "Moment", "moment must be positive");
    % Check maxItr
    maxItr = UseDefault(maxItr, maxItr >= 0, defaultMaxItr, "MaxItr", "maxItr must be nonnegative");
    % Check eigTol
    eigTol = UseDefault(eigTol, isUnitRange(eigTol), defaultEigTol, "EigTol", "eigTol must lie in [0,1)");
    % Check eigenspace
    if ~isempty(eigenspace) && (size(eigenspace, 1) ~= colA || ~all(isfinite(eigenspace), "all"))
        warning("CJSSRRR_hybrid:InvalidEigenspace", ...
                "eigenspace has size %d-by-%d; it must be finite with %d rows. Using default: [].", ...
                size(eigenspace,1), size(eigenspace,2), colA);
        eigenspace = defaultEigenspace;
    end
    % Check spectrum
    if ~isempty(spectrum)
        if numel(spectrum) < 2 || ~all(isfinite(spectrum)) || min(spectrum) == max(spectrum)
            warning("CJSSRRR_hybrid:InvalidSpectrum", ...
                    "spectrum must contain at least two finite values with distinct bounds. Using default: [].");
            spectrum = defaultSpectrum;
        else
            spectrum = [min(spectrum) max(spectrum)];
        end
    end
    % Check nev
    nev = UseDefault(nev, nev >= 0 && nev <= colA, defaultNev, "Nev", "nev must lie between 0 and size(A,2)");
    % Check V
    if ~isempty(V) && (size(V, 1) ~= colA || ~all(isfinite(V), "all"))
        warning("CJSSRRR_hybrid:InvalidV", ...
                "V has size %d-by-%d; it must be finite with %d rows. Using default: [].", ...
                size(V,1), size(V,2), colA);
        V = defaultV;
    end
    % Check degree
    degree = UseDefault(degree, degree >= 0, defaultDegree, "Degree", "degree must be nonnegative");
    % Check degreeD
    degreeD = UseDefault(degreeD, degreeD > 0, defaultDegreeD, "DegreeD", "degreeD must be positive");
    % Check omega
    omega = UseDefault(omega, omega > 0, defaultOmega, "Omega", "omega must be positive");
    if omega < 1
        warning("CJSSRRR_hybrid:TooSmallOmega", ...
                "Omega is less than 1, which may lead to failure in convergence.");
    end
    % Check restart
    restart = UseDefault(restart, any(strcmpi(restart, ["natural", "combine", "hybrid"])), defaultRestart, ...
                        "Restart", "restart must be 'natural', 'combine' or 'hybrid'");
    restart = lower(string(restart));
    % Check orthMethod
    orthMethod = UseDefault(orthMethod, any(strcmpi(orthMethod, ["TSVD", "QR"])), defaultOrthMethod, ...
                            "OrthMethod", "orthMethod must be 'TSVD' or 'QR'");
    orthMethod = upper(string(orthMethod));
    % Check svdTol
    svdTol = UseDefault(svdTol, isUnitRange(svdTol), defaultSvdTol, "SvdTol", "svdTol must lie in [0,1)");
    % Check eps_switch using the corrected eigTol
    defaultEps_switch = sqrt(eigTol);
    if any(strcmpi(parser.UsingDefaults, "eps_switch"))
        eps_switch = defaultEps_switch;
    else
        eps_switch = UseDefault(eps_switch, isUnitRange(eps_switch), defaultEps_switch, ...
                                "EpsSwitch", "eps_switch must lie in [0,1)");
    end
    % A single moment always uses QR for the effective orthonormalization
    if moment == 1
        orthMethod = "QR";
    end

end %~ end of CheckInputs


%   UseDefault: Warn and replace an invalid scalar or text option.
function value = UseDefault(value, isValid, defaultValue, parameter, reason)
    if isValid; return; end
    warning("CJSSRRR_hybrid:Invalid" + parameter, ...
            "%s. Using default: %s.", reason, string(defaultValue));
    value = defaultValue;
end %~ end of UseDefault



%   EstimateSpectrum: Estimate the spectrum of matrix A
function spectrum = EstimateSpectrum(A, step)
%   Ref: Y. Zhou and R.-C. Li, Bounding the spectrum of large Hermitian matrices,
%       Linear Algebra Appl., 435 (2011), pp. 480–493.
    % Lanczos procedure
    v = randn(size(A,2), 1); % Random initial vector
    [T, beta] = Lanczos(A, v, step);

    % Compute eigenvalues of the tridiagonal matrix T
    [V, d] = eig(T, "vector");
    d = real(d);

    % Estimate the spectrum
    l_min = d(1)   - beta * max(abs(V(end,1:2)));
    l_max = d(end) + beta * max(abs(V(end,end-1:end)));
    spectrum = [l_min l_max];
end %~ end of EstimateSpectrum


function [T, beta] = Lanczos(A, v, m)
%   Lanczos: m-step Lanczos procedure
    % Initialize
    T = zeros(m, m); % Tridiagonal matrix
    v = v / norm(v); % Normalize the initial vector

    % First step
    w      = A * v;
    T(1,1) = real(v'*w);     % First diagonal element
    w      = w - T(1,1) * v; % Orthogonalize
    beta   = norm(w);        % Compute beta

    % Loop for the remaining steps
    for i = 2:m
        T(i,i-1) = beta; % Sub-diagonal element
        T(i-1,i) = beta; % Super-diagonal element
        % Swap vectors
        tmp = w; w = v; v = tmp;
        % Current step
        v      = v / beta;         % Normalize
        w      = A * v - beta * w; % Modified Gram-Schmidt step
        T(i,i) = real(v'*w);       % Diagonal element
        w      = w - T(i,i) * v;   % Orthogonalize
        beta   = norm(w);
    end
end %~ end of Lanczos


%   EstimateDegree: Estimate the degree of the CJ polynomial series
function degree = EstimateDegree(a, b, moment, degree, degreeD)
%   Estimate the degree based on the interval and moments
%   degree = ceil( D * 2^(1/3) * pi^2 * moment^(1/6) / (b-a)^(4/3) ) - 2
    if degree > 0; return; end % If degree is already set, return
    degree = degreeD*2^(1/3)*pi^2*moment^(1/6) / (b-a)^(4/3);
    degree = ceil(degree) - 2;
end %~ end of EstimateDegree


%   GetCJCoeffs: Compute the CJ coefficients for the moment polynomials
function coeffs = GetCJCoeffs(a, b, degree, moment)
%   coeffs(j+1,k) = 2/π * \int_a^b cos(j*arccos(x)) * p(x, k) / \sqrt{1-x^2} dx
%                 = 2/π * \int_β^α cos(j*t) * p(cos(t), k) dt
%   The normalization factor 2/π will be replaced by 1/π for j = 0
    % Used polynomials
    T = @(t) cos((0:degree)' .* t);               % Trigonometric polynomials
    p = @(t, M) cos((M-1)*acos((2*t-a-b)/(b-a))); % Moment polynomial

    % Change the integration variable from x to arccos(x)
    a = acos(a); b = acos(b);
    % Initialize the coefficients matrix
    coeffs = zeros(degree+1, moment);
    % Moment 1 can be computed directly
    coeffs(1,1) = a-b;
    coeffs(2:end,1) = ( sin(a*(1:degree)) - sin(b*(1:degree)) )' ./ (1:degree)';
    % Compute the coefficients for higher moments by numerical integration
    for M = 2:moment
        coeffs(:,M) = integral(@(t) p(cos(t),M).*T(t), b, a, "ArrayValued", true);
    end

    % Multiply by the normalization factor
    coeffs(1,:)     = coeffs(1,:)     * (1/pi);
    coeffs(2:end,:) = coeffs(2:end,:) * (2/pi);

    % Jackson damping factors
    j = (0:degree)';
    theta = pi / (degree+2);
    damp = sin((j+1)*theta)/sin(theta) + (degree+1-j).*cos(j*theta);
    damp = damp / (degree+2);

    % Apply the damping factor to the coefficients
    coeffs = coeffs .* damp;
end %~ end of GetCJCoeffs


%   GetCJSubspace: Generate the CJ subspace based on the coefficients and starting matrix
function [S, vecTime, MVs_t, flops] = GetCJSubspace(A, c, r, V, coeffs, vecTime, MVs_t, flops)
    % V:       n-by-l
    % coeffs:  d-by-M
    % S:       n-by-(M*l)
    % vecTime: [MV-time flops-time]
    if isempty(V); S = []; return; end % If V is empty, return empty S

    % Count the MVs and other flops
    % MVs_t = d * l
    % flops = 2 * d * (1+M) * n * l
    degree = size(coeffs, 1) - 1;         % Degree of the CJ series
    moment = size(coeffs, 2);             % Number of moments
    MVs_t  = MVs_t + degree * size(V, 2); % Actual matrix-vector products
    flops  = flops + 2*degree*(1+moment)*numel(V);

    % Time recording
    time_MVs  = 0;
    time_flop = tic;

    % Map the spectrum interval of A to [-1, 1]
    A = (A - c*speye(size(A))) / r;

    % Generate the CJ matrix
    t  = tic;
    T0 = V; T1 = A * V;
    time_MVs = time_MVs + toc(t);

    % Scale A
    A = 2 * A;

    if moment == 1
        S = coeffs(1) * T0 + coeffs(2) * T1;
        for j = 3:degree+1
            % Compute the value of next Chebyshev polynomial
            t = tic;
            T = A * T1 - T0;
            time_MVs = time_MVs + toc(t);
            % Update T0, T1 and S
            T0 = T1; T1 = T;
            S  = S + coeffs(j) * T;
        end
    else % moment > 1
        S = coeffs(1,:) .* T0(:) + coeffs(2,:) .* T1(:);
        for j = 3:degree+1
            % Compute the value of next Chebyshev polynomial
            t = tic;
            T = A * T1 - T0;
            time_MVs = time_MVs + toc(t);
            % Update T0, T1 and S
            T0 = T1; T1 = T;
            S  = S + coeffs(j,:) .* T(:); % Faster than T(:)*coeffs(j,:)
        end
        % Reshape S to the correct size
        S = reshape(S, size(A, 1), []);
    end

    % Record the time
    time_flop = toc(time_flop) - time_MVs;
    vecTime = vecTime + [time_MVs time_flop];
end %~ end of GetCJSubspace


%   EstimateNev: Estimate the number of eigenvalues in the region
function [nev, S, V, timeOpts, MVs_t, flops] = EstimateNev(A, c, r, V, coeffs, timeOpts, MVs_t, flops)
    colV = size(V, 2); % Number of columns in V
    % Estimate the number of eigenvalues in the region
    time_nev = tic;
    [S, timeOpts.subspace, MVs_t, flops] = ...
        GetCJSubspace(A, c, r, V, coeffs, timeOpts.subspace, MVs_t, flops);
    timeOpts.nev = timeOpts.nev + toc(time_nev);
    nev = real(mean(sum(conj(V) .* S(:,1:colV), 1)));
end %~ end of EstimateNev


%   ResizeCJSubspace: Resize the CJ subspace to the desired size
function [S, V, timeOpts, MVs_t, flops] = ResizeCJSubspace(A, c, r, S, V, colV, coeffs, timeOpts, MVs_t, flops)
    colA    = size(A, 2);      % Number of columns in A
    colVold = size(V, 2);      % Number of columns in the old V
    moment  = size(coeffs, 2); % Number of moments
    % Resize V and S
    if colVold < colV
        U = randn(colA, colV - colVold); % Generate new random vectors
        [T, timeOpts.subspace, MVs_t, flops] = ...
            GetCJSubspace(A, c, r, U, coeffs, timeOpts.subspace, MVs_t, flops);
        V = [V U]; % Concatenate old and new V
        S = [S T]; % Concatenate old and new S
        ind_S = reshape((1:colVold)'      + colV*(0:moment-1), 1, []);
        ind_P = reshape((colVold+1:colV)' + colV*(0:moment-1), 1, []);
        S(:, [ind_S ind_P]) = S; % Reorder S to match the new V
    elseif colVold > colV
        % If V is larger than needed, truncate it
        V = V(:, 1:colV);
        ind = reshape((1:colV)' + colVold*(0:moment-1), 1, []);
        S = S(:, ind); % Reorder S to match the new V
    end
end %~ end of ResizeCJSubspace


%   RayleighRitz: Perform the standard Rayleigh-Ritz procedure
function [Y, E, AV, VAV] = RayleighRitz(A, V)
%   Standard Rayleigh-Ritz procedure for Hermitian matrices
    AV  = A * V;
    VAV = V' * AV;
    VAV = (VAV + VAV') / 2; % Ensure VAV is Hermitian
    % Compute the primitive Ritz pairs
    [Y, E] = eig(VAV, "vector");
    E = reshape(E, 1, []);  % Ensure E is a row vector
end %~ end of RayleighRitz


%   SinSpace: Compute sine between two subspaces
function angSin = SinSpace(V, U)
    if size(V,2) > size(U,2)
        s = size(V,2);
        [U, ~] = qr(U, "econ");
        W = [V U];
    else
        s = size(U,2);
        [V, ~] = qr(V, "econ");
        W = [U V];
    end
    R = qr(W, "econ");
    angSin = norm(R(s+1:end,s+1:end));
end %~ end of SinSpace


%  RefinedRitz: Compute refined Ritz pairs
function [Z, d, timeVec, r, idx_copy] = RefinedRitz(V, AV, VAV, rho, res, minGap, lout)

    t = tic;
    % Preprocess the parameters
    % The empty case is a valid iteration with no interval candidates
    if isempty(rho)
        Z = zeros(size(V,2), 0);
        d = zeros(1, 0);
        r = zeros(1, 0);
        idx_copy = false(1, 0);
        timeVec  = [toc(t) 0];
        return;
    end
    % Orthonormalize [V, AV] to accelerate the computation of SVD
    [H, J] = BGSOrth(V, AV, VAV);
    % Sort the residual norms
    [~, index] = sort(res, "ascend");
    rho = rho(index);


    % Initialize the parameters
    % The sizes of input matrices and vectors
    k = size(V, 2);     % The column size of V
    p = length(rho);    % The number of approximations
    % The primitive refined Ritz pairs
    Z = zeros(k, p);    % Primitive refined Ritz vectors
    d = zeros(1, p);    % Primitive refined Ritz values
    r = zeros(1, p);    % Residual of refined Ritz pairs
    % The SVD results
    Y = zeros(k, p, p); % The right singular vectors
    S = zeros(k, p);    % The singular values
    % The intermediate variable
    mult = 1:p;         % The i-th refined Ritz value has a multiplicity of mult(i)
    % The time record
    timeVec = [toc(t) 0];
    t = tic;


    % Compute the refined Ritz vectors
    for i = 1:p
        % Compute SVD
        [~, S(k:-1:1,i), Y(:,k:-1:1,i)] = svd(H - rho(i)*J, "econ", "vector");
        Z(:, i) = Y(:,1,i);
        d(   i) = real(Z(:,i)' * VAV * Z(:,i));
        r(   i) = norm(H*Z(:,i) - J*Z(:,i)*d(i));
    end


    % Reallocate the refined Ritz vectors
    idx_free = true(1, p);
    res_tmp  = max(r, minGap);
    for i = 1:p
        if ~idx_free(i); continue; end
        % Detect the multiplicity by residual norms
        idx_same = abs(d(i)-d) <= res_tmp(i) + res_tmp;
        idx_same = idx_same & idx_free;
        n_same   = nnz(idx_same);

        % Revise the multiplicity by singular values
        n_gap = nnz(S(1:n_same,i) > S(1,i)/eps^.5);
        if n_gap > 0
            idx_same(idx_same) = [true(1,n_same-n_gap) false(1,n_gap)];
            n_same = n_same - n_gap;
        end

        Z(:, idx_same) = Y(:,1:n_same,i);
        d(   idx_same) = real(sum(conj(Z(:,idx_same)).*(VAV*Z(:,idx_same)), 1));
        r(   idx_same) = vecnorm(H*Z(:,idx_same) - J*Z(:,idx_same).*d(idx_same));
        mult(idx_same) = i;
        idx_free = idx_free & (~idx_same);
    end


    % Remove extra/copy refined Ritz vectors by singular values
    idx_copy = false(1, p); % The indices of extra refined Ritz vectors
    res_tmp  = max(r, minGap);
    [~, ind] = sort(r, "descend");
    for i = ind
        idx_same = abs(d(i)-d) <= res_tmp(i) + res_tmp;

        ind_mult = unique(mult(idx_same));
        % If all the vectors are multiple, pass
        if length(ind_mult) <= 1; continue; end
        % Collect the multiple vectors
        for j = ind_mult
            idx_same = idx_same | mult==j;
        end

        % Detect extra refined Ritz vectors by singular values
        est_gap  = EstimateEigDis([d(~idx_same) lout], [r(~idx_same) 0 0], ...
                                d(idx_same&~idx_copy), r(idx_same&~idx_copy), minGap);
        idx_same = idx_same & ~idx_copy; % Remove the extra vectors
        if nnz(idx_same) == 0; continue; end
        inn_copy = GetCopy(Z(:,idx_same), d(idx_same), r(idx_same), est_gap);
        idx_copy(idx_same) = idx_copy(idx_same) | inn_copy;
    end

    timeVec(2) = toc(t);
end %~ end of RefinedRitz

%   BGSOrth: The block Gram-Schmidt orthonormalization
function [H, J] = BGSOrth(V, AV, VAV)
    [n, k] = size(V);

    % No orthonormalization when n is not large
    if n <= 2*k; H = AV; J = V; return; end

    % Block Gram-Schmidt orthonormalization
    H = [VAV; zeros(k,k)];
    W = AV - V * VAV;
    [W, H(k+1:end,:)] = qr(W, "econ");
    h = V' * W;
    W = W - V * h;
    R = qr(W, "econ");
    H(1:k,:) = H(1:k,:) + h * H(k+1:end,:);
    H(k+1:end,:) = R * H(k+1:end,:);
    J = speye(size(H));
end %~ end of BGSOrth

%   GetCopy: Detect extra refined Ritz vectors
function [idx_copy, U, S, V, delta] = GetCopy(Z, d, r, gap)
    idx_copy = false(size(d));
    eig_sin = r ./ gap;
    S = norm(Z(:,1));
    U = Z(:,1) / S;
    V = 1;
    delta = eig_sin(1)^2;

    for j = 2:length(d)
        [U1, S1, V1] = svdappend(U, S, V, Z(:,j));
        % Marked as extra if the singular value is too small
        if S1(end) <= max((delta + eig_sin(j)^2)^.5, S1(1)*128*eps)
            idx_copy(j) = true;
            continue;
        end
        % Update if the singular value is not small
        S = S1;
        U = U1;
        V = V1;
        delta = delta + eig_sin(j)^2;
    end
end %~ end of GetCopy

%   EstimateEigDis: Estimate distances between approximate values
function eigDis = EstimateEigDis(rho, res, d, r, minGap)
    % rho: the auxiliary approximate eigenvalues
    % res: the residual norms of rho
    % d:   the target approximate eigenvalues
    % r:   the residual norms of d

    if isempty(d); eigDis = []; return; end

    dis = abs(d - rho') - r - res.';
    dis(dis<=0) = Inf;
    dis = dis + r + res.';
    dis(dis<=minGap) = Inf;   % Set distances smaller than eigGap to Inf
    eigDis = min(dis, [], 1); % Minimum distance to other eigenvalues
end %~ end of EstimateEigDis
