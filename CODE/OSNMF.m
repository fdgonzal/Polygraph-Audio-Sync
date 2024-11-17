function [B, G, D] = OSNMF(X, K, alpha, beta, I)

% OSNMF - Performs Orthogonal and Sparse Non-Negative Matrix Factorization (OSNMF).
%
% Parameters:
%   X        - Input non-negative matrix (m x n).
%   K        - Number of components (rank of the factorization).
%   alpha    - Weight for the orthogonality penalty.
%   beta     - Weight for the sparsity penalty.
%   I        - Number of iterations for the optimization.
%
% Returns:
%   B - Basis matrix.
%   G - Activation matrix.
%   D - Divergence at each iteration.
%
% Reference:
%   Polygraph and Audio Synchronization Applied to Apnea Event Analysis 
%   Based on Non-negative Matrix Factorization (2024).
%

[m, n] = size(X);

% Initialize variables
D = zeros(1, I);  % Divergence values for each iteration
B = rand(m, K);   % Basis matrix
G = rand(K, n);   % Activation matrix
BG = B * G;       % Reconstruction

% Define penalty functions
phi_o = @(B) sum(sum(abs(B' * B - eye(K)))); % Orthogonality penalty
phi_s = @(G) sum(sum(abs(G * G' - eye(K))));      % Sparsity penalty

for iter = 1:I
    % Update basis matrix B
    num_B = (X ./ BG) * G';
    den_B = ones(m, n) * G' + 2 * alpha * phi_o(B) / (K * (K - 1));
    B = B .* (num_B ./ den_B);

    % Update activation matrix G
    den_G = BG + 2 * beta * phi_s(G) / (K * (K - 1));
    num_G = X ./ den_G;
    num_G = B' * num_G;
    G = sqrt(G .* num_G);

    % Update reconstruction
    BG = B * G;

    % Compute divergence
    D(iter) = sum(sum(X .* log(X ./ BG) - X + BG)) + alpha * phi_o(B) + beta * phi_s(G);
end

end
