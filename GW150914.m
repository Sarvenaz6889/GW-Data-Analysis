%  GW150914 — Time Series Plots of Fig. 1
%  Reproduces the three time-series rows of Fig. 1 in Abbott et al.,
%  PRL 116, 061102 (2016):
%     Row 1  H1 strain; L1 strain with H1 shifted and inverted on top
%     Row 2  the same data with the numerical relativity waveform over it
%     Row 3  residuals: data minus the filtered NR waveform

clear;
close all; 
clc;
% Initialization and Settings 
% Define paths for data and output figures.
dataDir = fullfile(getenv('HOME'),'Desktop','GW150914','data');
figDir  = fullfile(getenv('HOME'),'Desktop','GW150914','figs');
if ~exist(figDir,'dir'),
    mkdir(figDir); 
end
% Define input filenames (HDF5 format for raw strain, TXT for model)
fileH1 = fullfile(dataDir,'H-H1_LOSC_4_V2-1126259446-32.hdf5');
fileL1 = fullfile(dataDir,'L-L1_LOSC_4_V2-1126259446-32.hdf5');
fileNR = fullfile(dataDir,'fig2-unfiltered-template-reconstruction-H.txt');
% Long files, used only for the PSD estimate (see below)
fileH1long = fullfile(dataDir,'H-H1_LOSC_4_V2-1126257414-4096.hdf5');
fileL1long = fullfile(dataDir,'L-L1_LOSC_4_V2-1126257414-4096.hdf5');
% Physical and signal processing parameters
fs  = 4096;          % Sampling frequency in Hz
tZero  = 1126259462;    % GPS time of the figure's time origin
delay  = 0.0069;        % The GW reached Livingston 6.9 ms before Hanford
fLow  = 35;            % Bandpass lower cutoff frequency (Hz)
fHigh  = 350;           % Bandpass upper cutoff frequency (Hz)
% Tukey window roll-off, in seconds. The Guide specifies 0.5 s transition
% regions. Expressed in seconds rather than as a fraction so it does not
% change when the record length does.
rollSec = 0.5;
% PSD estimate span, in seconds. The Guide applies a Welch average to
% 1024 s of data centred on GPS 1126259462, in 4 s chunks spaced by 2 s.
psdSpan = 1024;
% Known instrumental noise lines to filter out (60 Hz power mains + harmonics)
% The 330 Hz is a known calibration line in the LIGO detectors.
notchHz = [60 120 180 240 300 330];
% Define the specific time window (in seconds) to plot, matching the
% published axes.
xLo = 0.285;
xHi = 0.46;
% Read the Strain Data 
% Extract the raw strain data array from the HDF5 files
h1 = h5read(fileH1,'/strain/Strain');
l1 = h5read(fileL1,'/strain/Strain');
% Read the starting GPS time from the file metadata.
% NOTE: Xstart is an attribute of /strain/Strain, not of /meta. Both
% groups exist in these files; reading the wrong one gives a wrong axis.
gps0 = double(h5readatt(fileH1,'/strain/Strain','Xstart'));
% Create a time vector relative to the figure's time origin. The caption
% gives times relative to September 14, 2015 09:50:45 UTC, which is GPS
% 1126259462 exactly — hence the published axis runs 0.30 to 0.45 s and
% the merger sits near 0.42 rather than at zero.
t = gps0 + (0:numel(h1)-1)'/fs - tZero;
% Load the theoretical waveform model provided by the LIGO collaboration.
% Column 1 is time in the same convention as above, column 2 is the NR
% strain scaled by 1e21.
NR = readmatrix(fileNR, 'CommentStyle','#');
tNR  = NR(:,1);             % Time vector for the model
hNR  = NR(:,2) * 1e-21;     % Scale up to physical strain units
% Apply a Tukey-like taper to the first 20ms of the theoretical waveform.
% The file covers only 0.25-0.45 s and is zero outside, so a hard step at
% 0.25 s would ring through the bandpass.
nTaper = round(0.020*fs);
ramp   = 0.5*(1 - cos(pi*(0:nTaper-1)'/nTaper));
taper  = ones(size(hNR));
taper(1:nTaper) = ramp;
hNR = hNR .* taper;
% Project the NR waveform onto the exact time grid of our raw data
nrH1 = interp1(tNR, hNR, t, 'linear', 0);
% Create the L1 theoretical waveform. Because L1 and H1 are oriented 
% differently and physically separated, we must shift the waveform 
% backward by the 6.9 ms travel time and invert it (multiply by -1).
% NOTE: GWOSC publishes the NR waveform projected onto H1 only. This is
% an approximation to the published L1 curve, which comes from the full
% projection onto L1's antenna pattern, so L1 agreement is correspondingly
% worse than H1's.
nrL1 = -interp1(tNR, hNR, t + delay, 'linear', 0);
% Noise Spectrum Estimation (PSD) 
% The Guide estimates the spectrum from 1024 s of data centred on the
% event, broken into overlapping 4 s chunks spaced by 2 s, each Tukey
% filtered. A longer average gives a smoother PSD than the 32 s file can.
% If the long files are absent, fall back to the 32 s record and say so.
if exist(fileH1long,'file')==2 && exist(fileL1long,'file')==2
    [psdH1, freq] = longPSD(fileH1long, tZero, psdSpan, fs);
    [psdL1, ~   ] = longPSD(fileL1long, tZero, psdSpan, fs);
    fprintf('PSD from %g s of data, as the Guide specifies.\n\n', psdSpan);
else
    [psdH1, freq] = pwelch(h1, hann(4*fs), 2*fs, 4*fs, fs);
    [psdL1, ~   ] = pwelch(l1, hann(4*fs), 2*fs, 4*fs, fs);
    warning(['Long files not found. Using the 32 s record for the PSD, ' ...
             'which is NOT what the Guide specifies.']);
end
% Whitening 
% Whiten the raw data AND the theoretical waveforms, identically.
% Whitening flattens the noise spectrum so low-frequency noise doesn't 
% artificially drown out the higher-frequency parts of the signal.
% Both must receive the same treatment, or the residual in row 3 would
% be dominated by the difference in processing rather than by any
% mismatch between model and data.
%
% NOTE: the Fig. 1 caption names only the bandpass and the notch filters.
% Following it literally, with no whitening, gives a figure with roughly
% twice the background of the published one, in which the signal is not
% clearly visible in the top row. The Guide's Section 4 supplies the
% missing step.
wh1 = whiten(h1,   freq, psdH1, fs, rollSec);
wl1 = whiten(l1,   freq, psdL1, fs, rollSec);
wnrH1 = whiten(nrH1, freq, psdH1, fs, rollSec);
wnrL1 = whiten(nrL1, freq, psdL1, fs, rollSec);
% Bandpass and Notch Filtering 
% Design an 8th-order IIR bandpass filter to isolate the 35-350 Hz band.
% designfilt's 'bandpassiir' default is Butterworth, matching the Guide.
% NOTE: the Guide points out that such narrow bandpassing is used only for
% visualization and is not employed in the LVC analyses.
bp = designfilt('bandpassiir','FilterOrder',8, ...
                'HalfPowerFrequency1',fLow, ...
                'HalfPowerFrequency2',fHigh, ...
                'SampleRate',fs);
% Apply the bandpass and notch filters. We use 'filtfilt' (zero-phase filtering)
% to run the filter forwards and backwards, ensuring the signal doesn't shift in time.
fh1   = applyFilters(wh1, bp, notchHz, fs);
fl1   = applyFilters(wl1, bp, notchHz, fs);
fnrH1 = applyFilters(wnrH1, bp, notchHz, fs);
fnrL1 = applyFilters(wnrL1, bp, notchHz, fs);
% Final Amplitude Scaling 
% Whitened data is in units of the noise standard deviation, not strain,
% so it cannot be labelled 1e-21 directly. Scale everything by one factor
% so the NR waveform peaks at 1.0 — the peak strain .
% The calibration comes from the model, not from the data, and the single
% factor leaves the relative amplitudes untouched.
scale = 1.0 / max(abs(fnrH1));
fh1   = fh1   * scale;
fl1   = fl1   * scale;
fnrH1 = fnrH1 * scale;
fnrL1 = fnrL1 * scale;
% Plotting the Three Time-Series Rows 
show = t > xLo-0.02 & t < xHi+0.02; % Define logical mask for the plot window
figure('Name','GW150914 Fig. 1 time series','Position',[80 80 1000 820]);
% ROW 1: Whitened, filtered strain data
subplot(3,2,1);
plot(t(show), fh1(show), 'r', 'LineWidth',1.1);
grid on;
xlim([xLo xHi]);
ylim([-1.5 1.5]);
ylabel('Strain (10^{-21})');
title('Hanford, Washington (H1)');
legend('H1 observed','Location','southwest');
subplot(3,2,2);
plot(t(show), fl1(show), 'b', 'LineWidth',1.1); 
hold on;
% Overlay H1 on L1 (shifted and inverted). The inversion accounts for the
% relative orientation of the two interferometers: the same passing wave
% stretches one arm pair where it squeezes the other.
% The shift is applied to the time axis, not with circshift, which would
% wrap the end of the array around to the start.
plot(t(show)-delay, -fh1(show), 'Color',[1 .5 .3], 'LineWidth',1.1);
grid on; 
xlim([xLo xHi]);
ylim([-1.5 1.5]);
title('Livingston, Louisiana (L1)');
legend('L1 observed','H1 observed (shifted, inverted)','Location','southwest');
% ROW 2: Data overlaid with the Numerical Relativity waveform
% NOTE: the published row 2 also shows shaded 90% credible regions for two
% independent waveform reconstructions. Those are credible bands from
% Bayesian inference, not available here, so only the NR waveform is drawn.
subplot(3,2,3);
plot(t(show), fh1(show), 'r', 'LineWidth',1.1); 
hold on;
plot(t(show), fnrH1(show), 'k', 'LineWidth',1.4);
grid on; 
xlim([xLo xHi]);
ylim([-1.5 1.5]);
ylabel('Strain (10^{-21})');
legend('H1 data','Numerical relativity','Location','southwest');
subplot(3,2,4);
plot(t(show), fl1(show), 'b', 'LineWidth',1.1);
hold on;
plot(t(show), fnrL1(show), 'k', 'LineWidth',1.4);
grid on;
xlim([xLo xHi]);
ylim([-1.5 1.5]);
legend('L1 data','Numerical relativity','Location','southwest');
% ROW 3: Residuals (Data MINUS Theoretical Waveform)
% If the model is a good description, the residual should look like noise.
resH1 = fh1 - fnrH1;
resL1 = fl1 - fnrL1;
subplot(3,2,5);
plot(t(show), resH1(show), 'r', 'LineWidth',1.0);
grid on;
xlim([xLo xHi]);
ylim([-0.6 0.6]);
xlabel('Time (s)');
ylabel('Strain (10^{-21})');
legend('Residual','Location','southwest');
subplot(3,2,6);
plot(t(show), resL1(show), 'b', 'LineWidth',1.0);
grid on;
xlim([xLo xHi]);
ylim([-0.6 0.6]);
xlabel('Time (s)');
legend('Residual','Location','southwest');
exportgraphics(gcf, fullfile(figDir,'fig1_timeseries.png'), 'Resolution',150);
% Residual RMS Test 
% NOT part of the published figure. Row 3 makes a claim by eye — that once
% the model is removed nothing structured is left. This is the same claim
% as a number: the residual in the signal window should fall to roughly
% the level of a stretch of pure noise.
inSignal = t > 0.35 & t < 0.45; % Time window containing the merger
inNoise  = t > 0.28 & t < 0.33; % Time window before the merger (pure noise)
fprintf(' RESIDUAL TEST \n\n');
fprintf('%-30s %8s %8s\n', '', 'H1', 'L1');
fprintf('%-30s %8.3f %8.3f\n', 'data rms, signal window', rms(fh1(inSignal)), rms(fl1(inSignal)));
fprintf('%-30s %8.3f %8.3f\n', 'residual rms, signal window', rms(resH1(inSignal)), rms(resL1(inSignal)));
fprintf('%-30s %8.3f %8.3f\n', 'residual rms, noise window', rms(resH1(inNoise)), rms(resL1(inNoise)));
fprintf('\n(units of 1e-21)\n\n');
%  Helper Functions
function [psd, freq] = longPSD(fname, tZero, span, fs)
% Welch PSD from a long stretch centred on the event, as the Guide
% specifies: 4 s chunks spaced by 2 s, each Tukey filtered. pwelch's
% hann window serves the same purpose here.
    strain = h5read(fname,'/strain/Strain');
    dt     = double(h5readatt(fname,'/strain/Strain','Xspacing'));
    gps0   = double(h5readatt(fname,'/strain/Strain','Xstart'));
    tRel   = gps0 + (0:numel(strain)-1)'*dt - tZero;
    sel    = tRel >= -span/2 & tRel <= span/2;
    x      = strain(sel);
    x      = x(~isnan(x));
    assert(numel(x) > 64*fs, 'PSD stretch too short in %s', fname);
    [psd, freq] = pwelch(x, hann(4*fs), 2*fs, 4*fs, fs);
end
function y = applyFilters(x, bp, notchHz, fs)
% Applies the bandpass filter, followed by notch filters for all specified 
% instrumental frequencies (like the 60Hz power grid noise).
    y = filtfilt(bp, x);
    for k = 1:numel(notchHz)
        f0 = notchHz(k);
        nf = designfilt('bandstopiir','FilterOrder',4, ...
                        'HalfPowerFrequency1',f0-1.5, ...
                        'HalfPowerFrequency2',f0+1.5, ...
                        'SampleRate',fs);
        y = filtfilt(nf, y); % filtfilt prevents phase shifting
    end
end
function w = whiten(x, freqAxis, psd, fs, rollSec)
% Divides the signal by the square root of the noise spectrum in the 
% frequency domain. This ensures that frequencies with high background noise 
% are suppressed, resulting in a flat noise spectrum with unit variance.
    x = x(:);
    n = numel(x);
    
    % Taper the edges of the data using a Tukey window to prevent the FFT 
    % from treating the start/end discontinuity as high-frequency noise.
    % rollSec is the transition region on EACH side, in seconds; tukeywin
    % takes the fraction of the whole record that is tapered, which is
    % twice that. The Guide specifies 0.5 s.
    r = min(2*rollSec*fs/n, 1);
    x = x .* tukeywin(n, r);
    X = fft(x);
    
    % Map the FFT bins to physical frequencies, folded at the Nyquist limit.
    f = (0:n-1)'*(fs/n);
    f = min(f, fs-f);
    
    % Interpolate the PSD to match the exact frequency bins of our FFT
    S = interp1(freqAxis, psd, f, 'linear', 'extrap');
    
    % Prevent division by zero or extremely small numbers. Strain PSD values 
    % are ~1e-46, so we must set a relative safety floor, not an absolute one.
    % max(S, eps) would replace EVERY value and silently disable whitening.
    S = max(S, min(psd(psd>0))*1e-6);
    
    % Divide the frequency-domain signal by sqrt(PSD) and transform back to time-domain.
    w = real(ifft(X ./ sqrt(S))) * sqrt(2/fs);
end