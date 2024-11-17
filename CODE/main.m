%% Automatic Synchronization Method Between Audio and Polygraph Signals for Sleep Studies
% 
% This script implements an automatic synchronization method to align audio 
% recordings with polygraph signals for sleep studies. Specifically, the
% method is based on the cross-correlation between an estimated respiratory
% sound signal anf the nasal flow signal from the polygraph.
%
% Reference:
%   Polygraph and Audio Synchronization Applied to Apnea Event Analysis 
%   Based on Non-negative Matrix Factorization (2024).
%
% -------------------------------------------------------------------------

%% Parameters

% Sampling frequencies
fs_audio = 4e3;  % Audio sampling rate (Hz).
fs_nf = 50;      % Nasal flow sampling rate (Hz).

% Filtering
cut_freq = 250;  % Cut-off frequency of the high-pass filter (Hz).

% STFT
N = 256;         % Window length (samples).
S = 0.5;         % Hop size (percentage of the window length).

% Segmentation
t_s = 180;       % Audio segment duration (seconds).
I_s = 4;         % Number of identical delays for algorithm termination.

% OSNMF
K = 20;          % Number of OSNMF components.
I = 300;         % Maximum OSNMF iterations.
alpha = 0.15;    % Orthogonality penalization weight.
beta = 0.1;      % Sparsity penalization weight.

% Delay estimation
t_max = 1200;    % Maximum delay between audio and polygraph recordings (seconds).
t_pow = 1;       % Time interval for power computation (seconds).

% File selection 
%[audio_file, audio_path] = uigetfile('*.wav', 'Select audio file');
%[nf_file, nf_path] = uigetfile('*.csv', 'Select nasal flow file');

t_sync = [];  % Obtained delay between audio and nasal flow signals.

%% Stage I: Preprocessing
audio_info = audioinfo([audio_path, audio_file]);
t_wav = audio_info.Duration;          % Total duration of the audio file (seconds).
fs_wav = audio_info.SampleRate;       % Sampling rate of the audio file (Hz).
J = fix(t_wav / t_s);                 % Number of segments in the audio file.
j = randperm(J);                      % Random permutation of segment indices.

time_delay = zeros(1, J);             % Array to store time delays for each segment.

for n_iter = 1:J  % Iterative algorithm
    %% Audio segment selection
    n1_wav = floor(fs_wav * (j(n_iter) - 1) * t_s + 1);
    n2_wav = floor(fs_wav * j(n_iter) * t_s);
    x_j = audioread([audio_path, audio_file], [n1_wav, n2_wav]);

    % Stereo to mono conversion
    if size(x_j, 2) == 2
        x_j = mean(x_j, 2);
    end

    % Resampling
    if fs_wav ~= fs_audio
        x_j = resample(x_j, fs_audio, fs_wav);
    end

    % High-pass filtering
    [b, a] = butter(6, cut_freq / (fs_audio / 2), 'high');
    x_j = filter(b, a, x_j);

    % STFT computation
    window = hamming(N);                     % Window function.
    nfft = 2^nextpow2(N * 2);                % Number of FFT points.
    hop_samples = round(S * N);             % Hop size in samples.
    noverlap = N - hop_samples;             % Overlap in samples.

    X = STFT(x_j, nfft, window, noverlap);   % Short-Time Fourier Transform.
    X_j = abs(X) / sum(sum(abs(X)));        % Normalized magnitude spectrogram.

    %% Stage II: OSNMF
    [B, G] = OSNMF(X_j, K, alpha, beta, I);

    %% Stage III: Optimal Component Selection
    SF = geomean(B.^2) ./ mean(B.^2);       % Spectral flatness.
    [~, k] = max(SF);                       % Optimal component index.
    B_opt = B(:, k);
    G_opt = G(k, :);
    X_opt = B_opt * G_opt;                  % Reconstructed spectrogram.

    % Estimated respiratory sound signal
    M_opt = X_opt ./ X_j;
    x_j_opt = iSTFT(M_opt .* X, nfft, window, noverlap);

    %% Stage IV: Delay Estimation

    % Read nasal flow interval
    x_p = readmatrix([nf_path, nf_file]);
    n1_nf = floor(fs_nf * j(n_iter) * t_s + 1);
    n2_nf = min(floor(fs_nf * ((j(n_iter) + 1) * t_s + t_max)), length(x_p));
    x_p = x_p(n1_nf:n2_nf);

    % Compute power signals
    P_xj = short_pow(x_j_opt, fs_audio, t_pow);
    P_xp = short_pow(x_p, fs_nf, t_pow);

    % Cross-correlation
    x_corr = xcorr1way_norm(P_xp, P_xj, t_max/t_pow + 1);
    [~, tau_max] = max(x_corr);            % Maximum cross-correlation index.
    time_delay(n_iter) = t_pow * (tau_max - 1); % Delay in seconds.

    % Check termination condition
    if sum(time_delay(n_iter) == time_delay(1:n_iter)) == I_s
        t_sync = time_delay(n_iter);       % Synchronized delay.
        break
    end

end
