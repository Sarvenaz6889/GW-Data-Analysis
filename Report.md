Project Report: Reproducing GW150914 and GW170817 Discoveries
Sarvenaz Sadeghhasani/ Date: October 2026
 
1. Introduction & Data
This report details the signal processing methods used to reproduce the time-series plots of GW150914 (Abbott et al., 2016, Fig. 1) and the spectrogram of GW170817 (Abbott et al., 2017, Fig. 2). Raw interferometer data and numerical relativity (NR) templates were sourced from GWOSC at a 4096 Hz sampling rate. The noise-subtracted release was used for GW170817.
 
2. General Signal Processing
To extract transient signals from dominant detector noise, the following steps were applied:
PSD Estimation. I estimated the noise floor using Welch's method in 4 s chunks spaced by 2 s. For GW150914 the average was taken over 1024 s centered on the event. For GW170817 I used a quiet window (−600 to −200 s) to prevent the 100 s signal from corrupting the noise estimate.
Whitening. I applied a Tukey window with a 0.5 s transition region on each side to prevent edge artifacts. Whitening was then achieved by dividing the Fourier coefficients by the ASD estimate to flatten the noise spectrum. The PSD must be floored before dividing, but the floor must be relative to the PSD: strain values are around 1e-46 while MATLAB's eps is 2.2e-16, so max(S, eps) replaces every value and silently disables the whitening without raising an error.
 
3. GW150914: Time-Domain Analysis
I applied a zero-phase, 8th order bandpass filter for the 35–350 Hz range. Narrow notch filters were added to remove 60 Hz power grid harmonics and the 330 Hz calibration line. The NR template was passed through an identical chain; otherwise, the residual would show the difference in processing rather than any model–data mismatch. The L1 signal was inverted and shifted by 6.9 ms to account for detector orientation and physical travel time. The GWOSC publishes the NR waveform projected onto H1 only, so L1 template is the H1 projection shifted and inverted.
 
4. GW170817: Time-Frequency Analysis
Spectrogram Generation. A Short-Time Fourier Transform (STFT) was computed using a 0.25 s Hann window with 95% overlap. The window length is a compromise: the early signal barely changes frequency and wants a long window, while the final second sweeps hundreds of Hz and wants a short one. No single STFT window handles both, which is why the published figures use a Q-transform.
Combination. I normalized each detector's spectrogram by its median to prevent the noisier detector from dominating, shifted L1 forward by 3 ms for phase alignment, and summed the arrays on a 20–500 Hz logarithmic scale. Tight color limits were needed: SNR 32.4 spread over 100 s is much fainter per pixel than GW150914's SNR 24 in 0.2 s.
 
References
1.	B. P. Abbott et al., Phys. Rev. Lett. 116, 061102 (2016).
2.	B. P. Abbott et al., Astrophys. J. Lett. 848, L12 (2017).
3.	B. P. Abbott et al., Class. Quantum Grav. 37, 055002 (2020), arXiv:1908.11170 — the Guide.
4.	Gravitational-Wave Open Science Center, https://gwosc.org

