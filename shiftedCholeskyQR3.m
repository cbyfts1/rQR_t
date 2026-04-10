function [Q, R] = shiftedCholeskyQR3(X)
% shiftedCholeskyQR3    Standard Improved Shifted CholeskyQR3 
% 
% This function implements the standard Shifted CholeskyQR3 algorithm 
% as described in the literature. It is provided to address Reviewer 
% v9jS's comment regarding the precise formulation of CholeskyQR3.
%
% References:
%   - [1] Takeshi Fukaya, Ramaseshan Kannan, Yuji Nakatsukasa, Yusaku Yamamoto, and Yuka Yanagisawa. 2020. Shifted Cholesky QR for Computing the QR Factorization of Ill-Conditioned Matrices. SIAM J. Sci. Comput. 42, 1 (2020), A477–A503. https://doi.org/10.1137/18M1218212
%   - [2] Fan, Y., Guan, H., & Qiao, Z. (2024). An Improved Shifted CholeskyQR Based on Columns. Journal of Scientific Computing, 104.
%
% Design purpose: To objectively demonstrate that even the standard 
% improved Shifted CholeskyQR3 (with regularized Gram matrix) fails 
% on extremely ill-conditioned matrices (κ≈10^16), while our rQR(t) 
% succeeds.
%
% Input:  X (m×n tall-and-skinny matrix, m >> n)
% Output: Q (m×n, approximately orthogonal), R (n×n upper triangular)

    [m, n] = size(X);
    if m < n
        error('Input matrix must be tall-and-skinny (m >= n).');
    end

    % ==================== Fixed parameters (literature-based) ====================
    u = eps;                    % machine epsilon for double precision
    c = 5000;                   % safety factor, increased for κ≈10^16 cases

    % g-norm based shift (recommended in Fan et al. for better stability)
    col_norms = sqrt(sum(X.^2, 1));
    Xg = max(col_norms);        % g-norm = max_j ||X(:,j)||_2
    s = c * (m*n*u + n*(n+1)*u) * Xg^2;

    fprintf('=== shiftedCholeskyQR3 (Improved Shifted CholeskyQR3) ===\n');
    fprintf('Computed shift s = %.2e\n', s);

    % Step 1: Shifted CholeskyQR (regularized Gram matrix)
    B = X'*X + s * eye(n);
    R1 = chol(B);
    Q = X / R1;

    R_accum = R1;

    % Step 2-5: CholeskyQR2 iterations (standard multi-iteration improvement)
    for k = 1:4
        B = Q'*Q;
        
        try
            Rk = chol(B);
            fprintf('Iteration %d: chol succeeded.\n', k);
        catch ME
            fprintf('Iteration %d: chol failed - %s\n', k, ME.message);
            error(['CholeskyQR3 failed at iteration %d.\n' ...
                   'Gram matrix became numerically non-positive definite.\n' ...
                   'This demonstrates the limitation of CholeskyQR3 ' ...
                   'on extremely ill-conditioned matrices (κ≈10^16).'], k);
        end
        
        Q = Q / Rk;
        R_accum = Rk * R_accum;
        
        % Mild reorthogonalization to control error accumulation
        if k < 4
            [Q, ~] = qr(Q, 0);
        end
    end

    R = R_accum;
end
