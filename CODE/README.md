# Algorithm for Automatic Synchronization of Audio and Polygraph Recordings in Sleep Studies

This repository implements an automatic synchronization method to align audio recordings with polygraph signals, specifically designed for sleep studies. The workflow includes preprocessing, source separation using Orthogonal Sparse Non-Negative Matrix Factorization (OSNMF), optimal component selection, and delay estimation through cross-correlation analysis. 

The repository provides MATLAB functions and scripts to process signals and calculate the synchronization delay, especially useful for a preliminary phase in the analysis of sleep disorders such as Obstructive Sleep Apnea (OSA).

## Files and Functionalities

### Core Functions

- **`STFT.m`**  
  Computes the Short-Time Fourier Transform (STFT) of a signal using a specified window, overlap, and FFT size. 

- **`iSTFT.m`**  
  Reconstructs a time-domain signal from its STFT matrix using the inverse Short-Time Fourier Transform (iSTFT).
  
- **`OSNMF.m`**  
  Implements the Orthogonal Sparse Non-Negative Matrix Factorization (OSNMF) algorithm to decompose a spectrogram into its components, including orthogonality and sparsity penalizations.

- **`short_pow.m`**  
  Computes the short-term power of a signal over specified time intervals, which is used to estimate respiratory activity in the nasal flow and audio signals.

- **`xcorr1way_norm.m`**  
  Calculates the one-way normalized cross-correlation between two signals, assuming one is delayed with respect to the other. Normalizes signals to have unit L2-norm before computation.

### Main Script

- **`synchronization_script.m`**  
  Main script to perform the synchronization process.
  
---
