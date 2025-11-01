function [P_solution, t_nodes] = MSA_RMDE_solver(A, B, C, R, F, t0, tf, N, epsilon, l)
% [P_solution, t_nodes] = MSA_RMDE_solver(A, B, C, R, F, t0, tf, N, epsilon, l)
%
% Solves the matrix Riccati differential equation (RMDE)
% with a terminal condition using the Multi-Step Simplified Algorithm (MSA, "Algorithm 4").
%
% The equation being solved is:
%    -dP/dt = A'P + PA - P(B*R_inv*B')P + C'C
%    P(tf) = F
%
% The function integrates backward from t_f to t_0.
%
% Input:
%   A         : (n x n) State matrix
%   B         : (n x s) Input matrix
%   C         : (p x n) Output matrix
%   R         : (s x s) Control weighting matrix (symmetric positive definite)
%   F         : (n x n) Terminal condition matrix (at t=tf)
%   t0, tf    : Integration start and end time (t0 < tf)
%   N         : Number of time steps
%   epsilon   : Precision tolerance for the algorithm
%   l         : Number of refinement iterations for the algebraic solution
%
% Output:
%   P_solution : (n x n x N) Solution matrix P(t) at N time nodes.
%                P_solution(:,:,j) is the solution at t_nodes(j).
%   t_nodes    : (1 x N) Time nodes corresponding to the solution (from t0 to tf - delta_T).
%
% Remarks:
%   This implementation is based on "Algorithm 4" and related equations (e.g., (23), (24)),
%   solving via backward recurrence.
%
% Example:
%   % (User should provide specific values for A, B, C, R, F, t0, tf, N, epsilon, l)
%   % [P, t] = MSA_RMDE_solver(A, B, C, R, F, 0, 10, 1000, 1e-9, 2);
%
%   Author                : [Ze Zhang]
%   E-mail                : [zezhang@smail.xtu.edu.cn]
%   Last modification     : [2025/11/01]

Q = C.' * C;
n_dim = size(A, 1);
I = eye(n_dim);

try
    R_inv = inv(R);
catch ME
    error('MSA_RMDE_solver:SingularR', 'Matrix R is singular. Cannot compute R_inv.');
end

BRB = B * R_inv * B.';
delta_T = (tf - t0) / N;

try
    iP_a_minus = care(-A, B, Q, R);
    P_a_minus = -iP_a_minus; % P_a_minus is the negative definite solution
catch ME
    warning('MSA_RMDE_solver:ARE_Failed', 'Failed to solve the Algebraic Riccati Equation (ARE).');
    rethrow(ME);
end

A_check = A - BRB * P_a_minus;

eig_A_check = eig(A_check);
abs_eig_A_check = abs(eig_A_check);
non_zero_eigs = abs_eig_A_check(abs_eig_A_check > eps);

if isempty(non_zero_eigs)
    warning("MSA_RMDE_solver:SingularA_check", "All eigenvalues of A_check are close to zero. Setting q = 1.0.");
    q = 1.0;
else
    max_eig_val = max(non_zero_eigs);
    min_eig_val = min(non_zero_eigs);
    q = sqrt(max_eig_val * min_eig_val);
end

qI = q * I;
try
    qI_plus_Acheck_inv = inv(qI + A_check);
catch ME
    fprintf("MSA_RMDE_solver:SingularMatrix", "(qI + A_check) matrix is singular (q = %f).\n", q);
    rethrow(ME);
end

A_hat = qI_plus_Acheck_inv * (qI - A_check);
A_hat_T = A_hat'; % Store A_hat^T
Q_hat = 2 * q * qI_plus_Acheck_inv * BRB * qI_plus_Acheck_inv';

rho_A_hat = max(abs(eig(A_hat)));
xi = (1 - rho_A_hat) / 2;
beta = (rho_A_hat + xi)^2;

norm_Q_hat_F = norm(Q_hat, 'fro');
if norm_Q_hat_F < eps
    norm_Q_hat_F = eps; % Avoid log(0)
end

log_num = log(epsilon * (1 - beta) / norm_Q_hat_F);
log_den = log(beta);
n_guess = floor(log_num / log_den) - l; % 'l' is the input parameter

if n_guess <= 0
    warning("MSA_RMDE_solver:N_Guess", "Calculated n_guess (%d) <= 0. Setting it to 1.", n_guess);
    n_guess = 1;
end

n = n_guess;
A_hat_n = A_hat^n;
condition_val = (rho_A_hat + xi)^n;

iter_count = 0;
max_iter = 10000; % Prevent infinite loop

while norm(A_hat_n, 'fro') >= condition_val
    n = n + 1;
    A_hat_n = A_hat_n * A_hat;
    condition_val = (rho_A_hat + xi)^n;

    iter_count = iter_count + 1;
    if iter_count > max_iter
        error("MSA_RMDE_solver:LoopError", "The while loop for determining 'n' exceeded %d iterations.", max_iter);
    end
end

m_0_prime = floor(sqrt(n)) + 1;
n_0_prime = floor(sqrt(n)) + 1;
% q_0_prime = 0; % As defined in the pseudo-code q_0' = 0

S_q0 = Q_hat;

S_n0_standard = compute_lyap_sum(A_hat, Q_hat, n_0_prime);
S_n0 = A_hat * S_n0_standard * A_hat_T;

A_n0_prime = A_hat^n_0_prime;
S_m0_standard = compute_lyap_sum(A_n0_prime, S_n0, m_0_prime-1);
S_m0 = A_n0_prime * S_m0_standard * A_n0_prime';

S_tilde_0 = S_q0 + S_n0 + S_m0;

P_hat_l = S_tilde_0;
for i = 1:l
    P_hat_l = A_hat * P_hat_l * A_hat_T + Q_hat;
end

t_nodes = t0 + (0:N-1) * delta_T;
P_solution = zeros(n_dim, n_dim, N);

E_term = P_hat_l;

X = expm(-A_check * delta_T);
X_T = X';

try
    P_next = inv(F - P_a_minus);
catch ME
    fprintf("MSA_RMDE_solver:SingularMatrix", "(F - P_a_minus) matrix is singular. Cannot compute recurrence initial value.\n");
    rethrow(ME);
end

for j = N:-1:1
    P_current = X * (P_next - E_term) * X_T + E_term;
    P_solution(:, :, j) = inv(P_current) + P_a_minus;
    P_next = P_current;
end

end

