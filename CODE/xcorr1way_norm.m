function x_corr = xcorr1way_norm(x, y, tau_max)

% xcorr1way_norm - Computes the one-way normalized cross-correlation between 
% two signals, considering only delays of y relative to x. Prior to cross-correlation,
% L2 normalization is applied.
%
% Parameters:
%   x        - Advanced signal (e.g., nasal flow) (vector).
%   y        - Delayed signal (e.g., audio) (vector).
%   tau_max  - Maximum shift (samples).
%
% Returns:
%   x_corr - Normalized cross-correlation between x and y (vector).
%
% Reference:
%   Polygraph and Audio Synchronization Applied to Apnea Event Analysis 
%   Based on Non-negative Matrix Factorization (2024).
%

% Normalize the delayed signal (y) to have unit L2-norm
y_norm = y / sqrt(sum(y.^2));

% Initialize cross-correlation output vector
x_corr = zeros(1, tau_max);

% Compute the normalized cross-correlation for each delay
for tau = 1:tau_max
    % Extract a segment of x corresponding to the current tau
    x_seg = x(tau:min(tau - 1 + length(y_norm), length(x)));
    x_seg_norm = x_seg / sqrt(sum(x_seg.^2));

    % Adjust lengths of the two signals to match
    min_len = min(length(x_seg_norm), length(y_norm));
    x_seg_norm = x_seg_norm(1:min_len);
    y_seg_norm = y_norm(1:min_len);

    % Compute the normalized cross-correlation
    x_corr(tau) = sum(x_seg_norm .* y_seg_norm) / min_len;
end

end

