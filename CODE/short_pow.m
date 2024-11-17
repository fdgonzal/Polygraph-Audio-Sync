function [x_pow] = short_pow(x, fs, t_pow)

% short_pow - Computes the short-term power of a signal over specified time intervals.
%
% Parameters:
%   x      - Input signal (vector).
%   fs     - Sampling frequency of the input signal (Hz).
%   t_pow  - Time interval for power calculation (seconds).
%
% Returns:
%   x_pow  - Short-term power of the signal (vector).
%
% Reference:
%   Polygraph and Audio Synchronization Applied to Apnea Event Analysis 
%   Based on Non-negative Matrix Factorization (2024).
%

% Determine the number of samples in an interval of t_pow seconds
inter_len = fs * t_pow; 
N = round(length(x) / inter_len); % Total number of intervals.

% Adjust interval sizes if they are not consistent
if rem(inter_len, 1) ~= 0
    inter_len   = floor(inter_len);
    dif_samples = length(x) - inter_len * N;
    inter_len   = [repmat(inter_len, 1, N - dif_samples), repmat(inter_len + 1, 1, dif_samples)];
end

% Preallocate the output vector for short-term power
x_pow = zeros(1, N);

% Compute the power for each interval
for j = 1:N  
    if isscalar(inter_len)
        x_frame = x(inter_len * (j - 1) + 1 : min(j * inter_len, length(x)));
    else
        x_frame = x(sum(inter_len(1:j-1)) + 1 : min(sum(inter_len(1:j)), length(x)));
    end
    % Calculate short-term power
    x_pow(j) = sum(x_frame .^ 2) / length(x_frame);
end

end
