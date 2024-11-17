function Y = STFT(x, nfft, window, noverlap)

% STFT - Computes the Short-Time Fourier Transform (STFT) of a signal.
%
% Parameters:
%   x        - Input signal (vector).
%   nfft     - Number of points for the FFT.
%   window   - Window function to be applied to each segment (vector).
%   noverlap - Number of overlapping samples between segments.
%
% Returns:
%   Y - STFT matrix containing complex STFT coefficients.
%
% Reference:
%   Polygraph and Audio Synchronization Applied to Apnea Event Analysis 
%   Based on Non-negative Matrix Factorization (2024).
%

nx = length(x);
nw = length(window);

% Zero-pad x if its length is less than the window length
if nx < nw
    x(nw) = 0;
    nx = nw;
end

x = x(:);
window = window(:);

% Calculate the number of segments for the STFT
ncols = fix((nx - noverlap) / (nw - noverlap));
cols_idx = 1 + (0:(ncols - 1)) * (nw - noverlap);
rows_idx = (1:nw)';

% Zero-pad x if necessary for the last segment
if length(x) < (nw + cols_idx(ncols) - 1)
    x(nw + cols_idx(ncols) - 1) = 0;
end

% Create the matrix of signal segments
seg = zeros(nw, ncols);
seg(:) = x(rows_idx(:, ones(1, ncols)) + cols_idx(ones(nw, 1), :) - 1);
seg = window(:, ones(1, ncols)) .* seg;

% Perform FFT on each segment
seg = fft(seg, nfft);

% Select frequency components
if isreal(x)
    if mod(nfft, 2) == 0
        sel = 1:(nfft / 2 + 1);
    else
        sel = 1:((nfft + 1) / 2);
    end
    Y = seg(sel, :);
else
    Y = seg;
end

end
