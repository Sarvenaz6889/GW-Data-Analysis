%  GW170817 — The Combined H1 + L1 Spectrogram
%  Reproduces the gravitational-wave inset of Fig. 2 in Abbott et al.,
%  ApJL 848, L12 (2017): the combined spectrogram of the signals received
%  by LIGO-Hanford (H1) and LIGO-Livingston (L1). Virgo is excluded
%  because its signal-to-noise ratio was too low to visibly contribute to
%  this specific visualization.
%
%  GW170817 (a binary neutron star merger) lasted for ~100 seconds in 
%  the detector sensitivity band.
%
%  Processing follows Section 4 of Abbott et al., "A guide to LIGO-Virgo
%  detector noise and extraction of transient gravitational-wave
%  signals", CQG 37, 055002 (2020) — the Guide.

clear;
close all; 
clc;

% Initial Settings & Parameters 
cfg.dataDir = '/Users/sarvenazsadeghhasani/Desktop/DATA';

% Filenames for the 4096 Hz, 2048-second datasets.
% CLN = the noise-subtracted data release from GWOSC.
cfg.fileH1 = 'H-H1_LOSC_CLN_4_V1-1187007040-2048.hdf5';
cfg.fileL1 = 'L-L1_LOSC_CLN_4_V1-1187007040-2048.hdf5';

cfg.fs     = 4096;
cfg.tEvent = 1187008882.4;    % GPS time of the merger

% Physical Arrival Delay: the signal reached Virgo first, Livingston
% 22 ms later, and Hanford 3 ms after that. Only the L1-to-H1 gap
% matters here, since Virgo is not in this figure.
cfg.dtLH = 0.003;

% Define the time window to plot relative to the merger. We only plot the 
% final 10 seconds because the early ~90 seconds are visually too faint.
cfg.tStart = -10;
cfg.tStop  =   2;

% Define the quiet window to calculate background noise (PSD).
% Must be far before the 100s signal to prevent the signal from whitening itself.
cfg.psdFrom = -600;
cfg.psdTo   = -200;

% Display frequency range
cfg.fMin = 20;
cfg.fMax = 500;

% STFT window length (seconds). 0.25s is a compromise: it tracks the fast 
% high-frequency sweep well, though it blurs the slow low-frequency start.
cfg.winSec = 0.25;

% Tukey window roll-off on EACH side, in seconds. The Guide specifies
% 0.5 s transition regions. Kept in seconds rather than as a fraction so
% it does not change when the segment length does.
cfg.rollSec = 0.5;

% Color limits to isolate the faint signal from the background noise
cfg.climSingle = [0 2.5];     % Limits for individual H1/L1 panels
cfg.climComb   = [0 4];       % Limits for the combined panel

% Read the Strain Data 
pathH1 = fullfile(cfg.dataDir, cfg.fileH1);
pathL1 = fullfile(cfg.dataDir, cfg.fileL1);

% Validate file existence to avoid obscure HDF5 internal errors later
assert(exist(pathH1,'file')==2, ...
    'H1 file not found: %s', pathH1);
assert(exist(pathL1,'file')==2, ...
    'L1 file not found: %s', pathL1);

[strainH1, tH1] = loadStrain(pathH1, cfg.tEvent);
[strainL1, tL1] = loadStrain(pathL1, cfg.tEvent);

fprintf('H1: %d samples, t = %.0f to %.0f s from merger\n', numel(strainH1), tH1(1), tH1(end));
fprintf('L1: %d samples, t = %.0f to %.0f s from merger\n\n', numel(strainL1), tL1(1), tL1(end));

% Noise Spectrum (PSD) Estimation
% Calculate the Power Spectral Density exclusively from the quiet window.
psdH1 = quietPSD(strainH1, tH1, cfg.psdFrom, cfg.psdTo, cfg.fs);
[psdL1, freq] = quietPSD(strainL1, tL1, cfg.psdFrom, cfg.psdTo, cfg.fs);

figure('Name','ASD','Position',[100 100 900 450]);
loglog(freq, sqrt(psdH1), 'r');
hold on;
loglog(freq, sqrt(psdL1), 'b');
grid on; 
xlim([10 2000]);
ylim([1e-24 1e-19]);
xlabel('Frequency (Hz)'); 
ylabel('ASD (strain/\surdHz)');
legend('H1','L1');
title('Noise spectrum, estimated prior to the signal');

% Sanity check on the downloads. H1 and L1 noise levels should be
% reasonably close in band for this event. If H1 sits far above L1,
% suspect a truncated download or a file that is not the CLN version.
bandOK = freq > 50 & freq < 300;
fprintf('Median ASD in 50-300 Hz:  H1 = %.2e   L1 = %.2e\n', ...
        median(sqrt(psdH1(bandOK))), median(sqrt(psdL1(bandOK))));
fprintf('(If H1 is much the larger, re-check the file sizes on disk.)\n\n');

% Segment Extraction and Whitening 
% Cut out the target analysis window with padding. Padding ensures that 
% the Tukey window (applied during whitening to prevent FFT ringing) 
% tapers the padding, not our physical signal.
pad = 8;    
[segH1, tSeg] = cutSegment(strainH1, tH1, cfg.tStart-pad, cfg.tStop+pad);
[segL1, ~   ] = cutSegment(strainL1, tL1, cfg.tStart-pad, cfg.tStop+pad);

% Whiten the data (divide by the square root of the PSD in frequency domain)
whH1 = whiten(segH1, freq, psdH1, cfg.fs, cfg.rollSec);
whL1 = whiten(segL1, freq, psdL1, cfg.fs, cfg.rollSec);

% Verify whitening success (Standard Deviation should be approximately 1)
fprintf('Whitened std:  H1=%.3f  L1=%.3f  (expect ~1)\n\n', std(whH1), std(whL1));

% Spectrogram Generation
winLen = round(cfg.winSec*cfg.fs);
nOver = round(0.95*winLen);   % 95% overlap ensures a visually continuous track
nFFT = 2^nextpow2(4*winLen); % Zero-padding for smoother visual interpolation

[Sh, f, tSpecH] = spectrogram(whH1, hann(winLen), nOver, nFFT, cfg.fs);
[Sl, ~, tSpecL] = spectrogram(whL1, hann(winLen), nOver, nFFT, cfg.fs);

% Shift the spectrogram time axes from absolute to merger-relative
tSpecH = tSpecH + tSeg(1);
tSpecL = tSpecL + tSeg(1);

% Physical Phase Alignment: Shift L1 forward by the 3 ms travel time.
% At this resolution the shift is below one pixel, but it is the
% physically correct thing to do and matters if the window is shortened.
tSpecL = tSpecL + cfg.dtLH;

keep = f >= cfg.fMin & f <= cfg.fMax;

% Median Normalization: Divide each array by its median noise level.
% This prevents the locally noisier detector from dominating the final sum.
% The median is used rather than the mean because most pixels are noise,
% so it measures the background and is barely moved by the signal.
Ah = abs(Sh(keep,:)); 
Ah = Ah ./ median(Ah(:));
Al = abs(Sl(keep,:));
Al = Al ./ median(Al(:));

% Interpolate L1 onto H1's exact time grid and sum them for the combined
% view. NOTE: "combined" is not defined in the Fig. 2 caption. This is
% one reading of it, and should be stated as a choice rather than as a
% reproduction of a specified method.
AlOnH = interp1(tSpecL, Al', tSpecH, 'linear', 0)';
Acomb = Ah + AlOnH;

% Plotting
figure('Name','GW170817 spectrogram','Position',[100 100 1000 800]);

% H1 Spectrogram
subplot(3,1,1);
imagesc(tSpecH, f(keep), Ah);
axis xy;
set(gca,'YScale','log'); % Log-scale 
ylim([cfg.fMin cfg.fMax]);
xlim([cfg.tStart cfg.tStop]);
clim(cfg.climSingle);
colormap(parula); 
colorbar;
ylabel('Frequency (Hz)');
title('H1 — Hanford');

% L1 Spectrogram
subplot(3,1,2);
imagesc(tSpecL, f(keep), Al);
axis xy; 
set(gca,'YScale','log');
ylim([cfg.fMin cfg.fMax]); 
xlim([cfg.tStart cfg.tStop]);
clim(cfg.climSingle); 
colormap(parula); 
colorbar;
ylabel('Frequency (Hz)');
title('L1 — Livingston');

% Combined Spectrogram
subplot(3,1,3)
imagesc(tSpecH, f(keep), Acomb);
axis xy;
set(gca,'YScale','log');
ylim([cfg.fMin cfg.fMax]); 
xlim([cfg.tStart cfg.tStop]);
clim(cfg.climComb);
colormap(parula); 
colorbar;
xlabel('Time from merger (s)'); 
ylabel('Frequency (Hz)');
title('Combined H1 + L1 — compare with the GW inset of Fig. 2');

figDir = fullfile(getenv('HOME'),'Desktop','GW170817','figs');
if ~exist(figDir,'dir'),
    mkdir(figDir); 
end
exportgraphics(gcf, fullfile(figDir,'fig2_spectrogram.png'), 'Resolution',150);

%  Helper Functions
function [strain, tRel] = loadStrain(fname, tEvent)
% Reads the HDF5 strain dataset and constructs a time axis centered on the merger.
% NOTE: Xspacing and Xstart are attributes of /strain/Strain, not of /meta.
% Both groups exist in these files; reading the wrong one gives a wrong axis.
    strain = h5read(fname,'/strain/Strain');
    dt     = double(h5readatt(fname,'/strain/Strain','Xspacing'));
    gps0   = double(h5readatt(fname,'/strain/Strain','Xstart'));
    n      = numel(strain);
    tRel   = gps0 + (0:n-1)'*dt - tEvent;
end

function [seg, tSeg] = cutSegment(strain, tRel, tFrom, tTo)
% Extracts the precise temporal segment required for analysis.
    sel  = tRel >= tFrom & tRel <= tTo;
    seg  = strain(sel);
    tSeg = tRel(sel);
    assert(~isempty(seg), 'Requested segment %g to %g s is outside the file.', tFrom, tTo);
end

function [psd, freq] = quietPSD(strain, tRel, tFrom, tTo, fs)
% Computes Welch's Power Spectral Density over a specified quiet window.
% NaNs are stripped rather than zeroed: a zero is a real value to the FFT
% and would ring across the spectrum.
    sel = tRel >= tFrom & tRel <= tTo;
    x   = strain(sel);
    x   = x(~isnan(x));
    assert(numel(x) > 8*fs, 'PSD stretch too short. Widen cfg.psdFrom / cfg.psdTo.');
    nfft = 4*fs;                      
    [psd, freq] = pwelch(x, hann(nfft), nfft/2, nfft, fs);
end

function w = whiten(x, freqAxis, psd, fs, rollSec)
% Whitens the signal by dividing the frequency-domain strain by the 
% square root of the PSD. This equalizes noise power across all frequencies.
    x = x(:);
    n = numel(x);
    
    % Tukey window gently tapers the boundaries to 0, preventing high-frequency 
    % spectral leakage when computing the FFT. rollSec is the transition
    % region on EACH side in seconds; tukeywin takes the fraction of the
    % whole record that is tapered, which is twice that.
    r = min(2*rollSec*fs/n, 1);
    x = x .* tukeywin(n, r);
    X = fft(x);
    
    % Construct the frequency bins to match the two-sided FFT output
    f = (0:n-1)'*(fs/n);
    f = min(f, fs-f);
    
    % Interpolate the background PSD onto the precise FFT frequency grid
    S = interp1(freqAxis, psd, f, 'linear', 'extrap');
    
    % Enforce a noise floor limit relative to the PSD.
    % Using pure `eps` here would disable whitening because physical strain PSD 
    % values (~1e-46) are vastly smaller than standard floating-point epsilon.
    S = max(S, min(psd(psd>0))*1e-6);
    
    % Normalize, drop the imaginary numeric noise, and scale for unit variance
    w = real(ifft(X ./ sqrt(S))) * sqrt(2/fs);
end