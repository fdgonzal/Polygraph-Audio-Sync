function x = iSTFT(X, nfft, window, noverlap)

% iSTFT - Reconstructs a time-domain signal from its STFT using the inverse Short-Time Fourier Transform (ISTFT).
%
% Parameters:
%   X        - Input STFT matrix (complex).
%   nfft     - Number of FFT points used in the STFT.
%   window   - Window function applied to each segment (vector).
%   noverlap - Number of overlapping samples between segments.
%
% Returns:
%   x - Reconstructed time-domain signal (vector).
%
% Reference:
%   Polygraph and Audio Synchronization Applied to Apnea Event Analysis 
%   Based on Non-negative Matrix Factorization (2024).
%

% Get the size of the STFT matrix: number of frequencies (F) and time frames (T)
[num_freqs, num_frames] = size(X);

% Extend the spectrogram to form the full symmetric spectrum
num_extra = nfft - num_freqs;
mirror_part = conj(X(num_extra + 1:-1:2, :));

% Construct the full spectrum for IFFT and perform IFFT for each frame
full_spectrum = [X; mirror_part];
time_frames = real(ifft(full_spectrum));

% Set the window length
window_len = length(window);
time_frames = time_frames(1:window_len, :);

% Configure overlapping windows
win_matrix = repmat(window, 1, num_frames);
step_size = window_len - noverlap;

% Initialize output signal and window sum vector for normalization
reconstructed_signal = zeros(num_frames * step_size + noverlap, 1);
window_sum = reconstructed_signal;

% Overlap-add reconstruction
for idx = 1:num_frames
    % Index range in output signal for current window
    range = ((idx - 1) * step_size + (1:window_len))';

    % Add current IFFT frame to the output signal
    reconstructed_signal(range) = reconstructed_signal(range) + time_frames(:, idx);

    % Accumulate window sum for normalization
    window_sum(range) = window_sum(range) + win_matrix(:, idx);
end

window_sum(1:round(window_len / 2)) = max(1, window_sum(1:round(window_len / 2)));
window_sum(end-round(window_len / 2):end) = max(1, window_sum(end-round(window_len / 2):end));

x = reconstructed_signal ./ window_sum;

end
