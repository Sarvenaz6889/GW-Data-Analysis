Reproducing GW150914 and GW170817
Sarvenaz Sadeghhasani | October 2026
Reproduces the time-series rows of Fig. 1 of GW150914 (Abbott et al., 2016) and the spectrogram inset of Fig. 2 of GW170817 (Abbott et al., 2017). MATLAB R2024b.
 
1. Data
From GWOSC at 4096 Hz. For GW150914 the 4096 s files are needed rather than the 32 s ones, since the conditioning chain requires 264 s of PSD training data. For GW170817 the noise-subtracted (CLN) release was used.
Both data quality masks are checked first. Note the inverted convention on the injection mask: the bit names read NO_CW_HW_INJ, so a cleared bit means that injection is present.
 
2. Data conditioning
Implements Prof. Mohanty's LALdatacondchain.m (IFODATACOND, bandId branch), with every parameter taken from LALdatacondparams_4kHz.yml.
Detrending is a FIR highpass designed with fir1, cutoff 25 Hz, order 8192, applied with fftfilt. FIR rather than IIR so the phase is linear and the delay is exactly half the order.
PSD estimation uses pwelch_median over 32 segments of 16 s spaced by 8 s, each weighted by a Tukey window with cosine fraction 0.5. The median rather than the mean, because a segment holding a glitch pulls the mean up but barely moves the median.
Whitening is a FIR filter designed with fir2 whose magnitude response is 1/√PSD, with a 16 s impulse response. The key step is that the PSD is levelled outside 25 to 1650 Hz rather than zeroed — held at the band-edge value. Below 25 Hz the data is already gone, so 1/√PSD would explode; above 1650 Hz there is nothing worth whitening.
Downsampling is resample by a factor of 2 with a 64th-order anti-alias filter, then a scale of √(2/f_s). Finally the chain's windowing step zeroes the first 18 s, which it treats as contaminated by filter settling.
The chain records its accumulated filter delay as fltrTShft rather than correcting the time axis, which suits a matched filter but not a figure with a published time axis. I therefore remove the true group delay, 9.000 s.
 
3. GW150914
The Fig. 1 caption names only a 35–350 Hz bandpass and band-reject filters. Whitening is not mentioned, but without it the background is roughly twice the published level and the signal is not clearly visible in the top row. I found this empirically; Sec. 4 of the guide (CQG 37, 055002) then confirmed whitening is part of the chain.
The NR template goes through an identical filter chain as the data; otherwise the residual would show the difference in processing rather than any model–data mismatch. L1 is inverted and shifted by 6.9 ms for detector orientation and travel time, applied to the time axis rather than with circshift, which would wrap the array. GWOSC publishes the NR waveform projected onto H1 only, so the L1 template is that projection shifted and inverted, and L1 agreement is correspondingly worse.
Row 3 of the published figure quantifying it, in units of 1e-21: the data rms in the signal window is 0.320 for H1 and 0.228 for L1, falling to 0.122 and 0.143 once the model is subtracted, against 0.126 and 0.100 in a pure-noise window. For H1 the residual matches noise to within 3 percent.
The conditioned standard deviation comes out at 1.207 for H1 and 1.210 for L1, not 1. The reason is that pwelch_median returns the raw median across segments with no bias correction.
 
4. GW170817
A Q-transform, as the published figure uses, plus an STFT for comparison. The Q-transform holds Q = f/Δf constant so its window shortens as frequency rises, which is what a chirp needs: GW170817 barely moves in frequency for most of its length, then sweeps hundreds of Hz in the final second, and no single STFT window serves both ends.
"Combined" is not defined in the caption. I normalize each detector's map by its own median so the noisier one does not dominate, shift L1 forward by 3 ms for arrival alignment, and sum on a 20–500 Hz logarithmic scale. Tight colour limits were needed: network SNR 32.4 spread over 100 s is much fainter per tile than GW150914's 24 in 0.2 s.
The full-band conditioned standard deviation came out at 1.378 for H1 and 2.962 for L1, with no drift across the record, so the cause is a constant scale rather than non-stationary noise. Restricted to the band the figures actually show, 30 to 400 Hz, the two agree: 0.703 and 0.707, a ratio of 1.01 against 2.15 across the full band. The excess lies above 1650 Hz, where the chain levels the PSD instead of whitening, so narrow lines in the kHz region survive unflattened. They have no effect on the displayed figures.
 
References
1.	Abbott et al., Phys. Rev. Lett. 116, 061102 (2016).
2.	Abbott et al., Astrophys. J. Lett. 848, L12 (2017).
3.	Abbott et al., Class. Quantum Grav. 37, 055002 (2020), arXiv:1908.11170.
4.	Chatterji et al., Class. Quantum Grav. 21, S1809 (2004).
5.	S. D. Mohanty, IFODATACOND, bandId branch.
6.	GWOSC, https://gwosc.org

