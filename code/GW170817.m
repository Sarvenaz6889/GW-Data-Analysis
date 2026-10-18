%  GW170817 — Fig. 2 combined spectrogram, using LALdatacondchain
clear;
close all;
clc;
dataDir = '/Users/sarvenazsadeghhasani/Desktop/DATA';

% CLN = the noise-subtracted data release from GWOSC.
fileH1 = fullfile(dataDir,'H-H1_LOSC_CLN_4_V1-1187007040-2048.hdf5');
fileL1 = fullfile(dataDir,'L-L1_LOSC_CLN_4_V1-1187007040-2048.hdf5');

figDir = fullfile(getenv('HOME'),'Desktop','GW170817','figs');
if ~exist(figDir,'dir');
    mkdir(figDir); 
end

% Event parameters 
fs  = 4096;
tEvent = 1187008882.4;       % GPS merger time

% Arrival delay: Virgo first, Livingston 22 ms later, Hanford 3 ms after
% that. Only the L1-to-H1 gap matters here.
dtLH = 0.003;

% Stretch to condition, relative to the merger. Must hold 264 s of
% training plus 18 s settling, ending well before the 100 s signal.
%
% anaStop must leave room at the END too. The last settleT = 18 s of the
% conditioned record is discarded, and the group delay shifts the axis
% back another 9 s, so the usable data ends at roughly anaStop - 18.
% With anaStop = 20 the valid span stopped at -7 s and the merger fell
% outside it, which is what left the right half of the figures blank.
anaStart = -700;
anaStop  =   50;

% LALdatacondparams_4kHz.yml
P.hiPassFltrCut     = 25;      % Hz
P.hiPassFltrOrdr    = 8192;
P.psdSegLen         = 16;      % s
P.psdSegOvrlp       = 8;       % s
P.psdNSeg           = 32;
P.psdSegWinParam    = 0.5;     % Tukey cosine fraction
P.whtnHiFreq        = 1650;    % Hz
P.whtFltrImpRspLen  = 16;      % s
P.dnSmplFac         = 2;
P.dnSmplFltrOrdr    = 64;
trnSegStrt = 0;                % s, as in test_LALdatacondchain_testspec.yml

% Conditioning switches 
dtrndON  = 1;
whtnON   = 1;
dnSmplON = 1;

% Display window 
dsp.tStart = -10;              % the early ~90 s are too faint to see
dsp.tStop  =   2;
dsp.fMin   =  20;
dsp.fMax   = 500;

% Q-transform 
% Q = 12 keeps
% enough frequency resolution for the long narrowband inspiral while
% letting the final sweep stay together.
qt.Q     = 12;
qt.nFreq = 150;

% STFT, for comparison 
st.winSec = 0.25;

% Colour limits 
% Lowered from [0 6] and [0 10]. The chirp here is faint per tile —
% network SNR 32.4 spread over 100 s is much weaker per pixel than
% GW150914's SNR 24 in 0.2 s — so the ceiling has to sit close to the
% track's own amplitude or it disappears into the background.
climSingle = [0 3];
climComb   = [0 5];

%  Load and validate
[rawH1, tRawH1] = loadStrain(fileH1, tEvent, fs);
[rawL1, tRawL1] = loadStrain(fileL1, tEvent, fs);

fprintf('H1: %d samples, t = %.0f to %.0f s from merger\n', ...
        numel(rawH1), tRawH1(1), tRawH1(end));
fprintf('L1: %d samples, t = %.0f to %.0f s from merger\n\n', ...
        numel(rawL1), tRawL1(1), tRawL1(end));

validateDataQuality(fileH1, rawH1, 'H1', tRawH1);
validateDataQuality(fileL1, rawL1, 'L1', tRawL1);

selH1 = tRawH1 >= anaStart & tRawH1 <= anaStop;
selL1 = tRawL1 >= anaStart & tRawL1 <= anaStop;
datH1 = rawH1(selH1);
datL1 = rawL1(selL1);
tAna0 = tRawH1(find(selH1,1));

% The two releases may be centred differently, so trim to a common length.
nAna  = min(numel(datH1), numel(datL1));
datH1 = datH1(1:nAna);
datL1 = datL1(1:nAna);
clear rawH1 rawL1 tRawH1 tRawL1

tAna = tAna0 + (0:nAna-1)'/fs;
fprintf('Conditioning %.0f s, t = %+.0f to %+.0f s\n\n', nAna/fs, tAna(1), tAna(end));


%  STEP 1 — DETRENDING
if dtrndON
    % fir1 highpass applied causally with fftfilt. FIR rather than IIR,
    % so the phase is linear and the delay is exactly half the order.
    hpB = fir1(P.hiPassFltrOrdr, P.hiPassFltrCut/(fs/2), 'high');

    tic;
    dtH1 = fftfilt(hpB, datH1);
    dtL1 = fftfilt(hpB, datL1);
    fprintf('Step 1  detrending: fir1 highpass %g Hz, order %d  (%.1f s)\n', ...
            P.hiPassFltrCut, P.hiPassFltrOrdr, toc);
    hpDly = P.hiPassFltrOrdr/fs;
else
    dtH1 = datH1; 
    dtL1 = datL1; 
    hpDly = 0;
    fprintf('Step 1  detrending: OFF\n');
end

%  STEP 2 — PSD ESTIMATION
% Median-averaged Welch on a training segment. The median rather than the
% mean: each periodogram bin is exponentially distributed about the true
% PSD, so a segment containing a glitch pulls the mean up but barely
% moves the median.
nSegSmpl   = floor(P.psdSegLen   * fs);
nOvrlpSmpl = floor(P.psdSegOvrlp * fs);

strtTrnIdx = max(1, floor((trnSegStrt + hpDly)*fs));
% LALdatacondchain computes this end index from the start of the
% RECORD, not from the start of the training segment:
%     endTrnPSDIndx = min([nSmpls, (nSeg-nOvrlp)*(psdNSeg-1)+nSeg])
% so when trnSegStrt or the filter delay pushes the start forward,
% the training span is shorter than psdNSeg full segments. That is
% reproduced here rather than corrected.
endTrnIdx  = min(nAna, (nSegSmpl - nOvrlpSmpl)*(P.psdNSeg-1) + nSegSmpl);
assert(endTrnIdx > strtTrnIdx + 3*nSegSmpl, ...
       'Training span too short for a median PSD. Check trnSegStrt.');

trnSpan = [tAna(strtTrnIdx), tAna(endTrnIdx)];
trnLen  = diff(trnSpan);
needLen = ((nSegSmpl - nOvrlpSmpl)*(P.psdNSeg-1) + nSegSmpl)/fs;

fprintf('Step 2  PSD training: t = %+.1f to %+.1f s  (%.0f s, need %.0f s)\n', ...
        trnSpan(1), trnSpan(2), trnLen, needLen);

assert(trnLen >= 0.99*needLen, 'Training segment truncated. Move anaStart earlier.');

% CRITICAL for this event. The signal occupies roughly the last 100 s
% before the merger, so the training window must end well before -100 s.
assert(trnSpan(2) < -150, ...
       ['Training segment reaches the 100 s signal. Move anaStart ' ...
        'earlier or reduce psdNSeg.']);
fprintf('        margin to the signal: %.0f s\n', -100 - trnSpan(2));
psdWin = tukeywin(nSegSmpl, P.psdSegWinParam);
tic;
[pxxH1, pxxF] = pwelchMedian(dtH1(strtTrnIdx:endTrnIdx), nSegSmpl, nOvrlpSmpl, psdWin, fs);
[pxxL1, ~   ] = pwelchMedian(dtL1(strtTrnIdx:endTrnIdx), nSegSmpl, nOvrlpSmpl, psdWin, fs);
fprintf('        %d segments of %g s, %g s overlap, Tukey r = %g  (%.1f s)\n', ...
        P.psdNSeg, P.psdSegLen, P.psdSegOvrlp, P.psdSegWinParam, toc);
fprintf('        frequency resolution %.4f Hz\n', 1/P.psdSegLen);

figure('Name','ASD','Position',[100 100 900 450]);
loglog(pxxF, sqrt(pxxH1), 'r');
hold on
loglog(pxxF, sqrt(pxxL1), 'b');
grid on;
xlim([10 2000]);
ylim([1e-24 1e-19]);
xlabel('Frequency (Hz)');
ylabel('ASD (strain/\surdHz)');
legend('H1','L1');
title('Median-averaged ASD, from before the signal');

% Sanity check on the downloads. H1 and L1 should be reasonably close in
% band for this event. If H1 sits far above L1, suspect a truncated
% download or a file that is not the CLN version.
bandOK = pxxF > 50 & pxxF < 300;
fprintf('        median ASD 50-300 Hz:  H1 = %.2e   L1 = %.2e\n', ...
        median(sqrt(pxxH1(bandOK))), median(sqrt(pxxL1(bandOK))));

%  STEP 3 — WHITENING
% The PSD is LEVELLED outside [hiPassFltrCut, whtnHiFreq]: held at the
% band-edge value rather than used as measured. Below 25 Hz the data is
% already gone, so 1/sqrt(PSD) would explode; above 1650 Hz there is
% nothing worth whitening. Levelling gives a flat finite response in both
% regions and keeps the fir2 design well conditioned.
whtFltrOrdr = floor(P.whtFltrImpRspLen * fs);

if whtnON
    fprintf('Step 3  whitening: fir2, %g s impulse response (order %d)\n', ...
            P.whtFltrImpRspLen, whtFltrOrdr);
    fprintf('        PSD levelled below %g Hz and above %g Hz\n', ...
            P.hiPassFltrCut, P.whtnHiFreq);
    fprintf('        designing and applying (this is the slow step)...\n');
    tic;
    whtB_H1 = designWhitener(pxxH1, pxxF, fs, P, whtFltrOrdr);
    whtB_L1 = designWhitener(pxxL1, pxxF, fs, P, whtFltrOrdr);

    wH1 = fftfilt(whtB_H1, dtH1);
    wL1 = fftfilt(whtB_L1, dtL1);
    fprintf('        done (%.1f s)\n', toc);
else
    wH1 = dtH1;  wL1 = dtL1;
    fprintf('Step 3  whitening: OFF\n');
end

%  STEP 4 — DOWNSAMPLE
if dnSmplON
    fsD = fs / P.dnSmplFac;
    cH1 = resample(wH1, 1, P.dnSmplFac, P.dnSmplFltrOrdr);
    cL1 = resample(wL1, 1, P.dnSmplFac, P.dnSmplFltrOrdr);

    % LALdatacondchain's normalization: makes data at the noise floor
    % have unit variance.
    nrmFac = sqrt(2/fsD);
    cH1 = nrmFac * cH1;
    cL1 = nrmFac * cL1;

    fprintf('Step 4  downsampled %g -> %g Hz (filter order %d), scaled by %.4g\n\n', ...
            fs, fsD, P.dnSmplFltrOrdr, nrmFac);
else
    fsD = fs;  nrmFac = sqrt(2/fsD);
    cH1 = nrmFac*wH1;  cL1 = nrmFac*wL1;
    fprintf('Step 4  downsampling: OFF\n\n');
end


%Time axis, with the group delay removed 
% fftfilt is causal, so both FIR filters delay the output. For a
% linear-phase FIR of order N the group delay is exactly N/2 samples.
% resample() compensates its own delay internally.
grpDly  = (P.hiPassFltrOrdr*dtrndON + whtFltrOrdr*whtnON) / (2*fs);
settleT = (P.hiPassFltrOrdr*dtrndON + whtFltrOrdr*whtnON) / fs;

% Data windowing, as LALdatacondchain does
% The chain zeroes everything before the accumulated filter delay rather
% than trimming it, and records that delay as fltrTShft in the output
% file so a downstream consumer knows the data is shifted:
strtTCutIdx = min(numel(cH1), floor(max(0,settleT)*fsD) + 1);
nZeroed     = strtTCutIdx - 1;
cH1(1:nZeroed) = 0;
cL1(1:nZeroed) = 0;

t = tAna0 + (0:numel(cH1)-1)'/fsD - grpDly;
fprintf('Group delay removed: %.3f s   (settling time %.3f s)\n', grpDly, settleT);
valid = t > (tAna0 + settleT - grpDly) & t < (t(end) - settleT);
fprintf('Zeroed the first %d samples (%.1f s), as the chain does.\n', ...
        nZeroed, nZeroed/fsD);
fprintf('Valid span: t = %+.1f to %+.1f s\n\n', ...
        t(find(valid,1)), t(find(valid,1,'last')));

% Conditioning sanity check 
% Conditioned data should have standard deviation near 1 away from the
% signal. Far from 1 means the PSD does not describe this stretch — most
% often because the training window overlapped a period of different
% detector behaviour, which the DQ report above would show.
quiet  = valid & t < -200;
stdH1  = std(cH1(quiet));
stdL1  = std(cL1(quiet));
fprintf('Conditioned std (quiet region):  H1 = %.3f   L1 = %.3f   (expect ~1)\n', ...
        stdH1, stdL1);

for s = [stdH1 stdL1]
    if s < 0.8 || s > 1.25
        warning(['Conditioned std is %.2f, not ~1. See the block scan ' ...
                 'below: if the value drifts across the record the noise ' ...
                 'is non-stationary and the training window is simply too ' ...
                 'far from the signal.'], s);
    end
end

% Block-by-block variance scan.
%
% The DQ report has already ruled out the obvious explanation. for this
% event the only flagged seconds sit on the merger itself. 
% The signature to look for: a value near 1 at the start of the record,
% drifting away as t increases. That says move the training window
% closer to the signal.
fprintf('\nWhitened std in 50 s blocks (watch for drift):\n');
edges = -650:50:-100;
for k = 1:numel(edges)-1
    blk = valid & t > edges(k) & t < edges(k+1);
    if any(blk)
        fprintf('  %+5.0f to %+5.0f s:   H1 = %.2f   L1 = %.2f\n', ...
                edges(k), edges(k+1), std(cH1(blk)), std(cL1(blk)));
    end
end
fprintf('\n');

% In-band variance check 
% The standard deviation above is measured across the FULL band, up to
% Nyquist. But whtnHiFreq levels the PSD above 1650 Hz instead of
% whitening there, so any narrow lines in that region survive unflattened
% and inflate the total variance — without affecting the 20-500 Hz region
% the spectrograms display.
%
% For this event the printed median ASD showed L1 QUIETER than H1 in
% 50-300 Hz, yet its whitened variance came out roughly twice as large.
% That combination points at out-of-band power rather than a problem with
% the conditioning. Measuring inside the display band settles it: if the
% two detectors agree here, the full-band figures are a red herring.
bpChk = designfilt('bandpassiir','FilterOrder',8, ...
                   'HalfPowerFrequency1',30, ...
                   'HalfPowerFrequency2',400, ...
                   'SampleRate',fsD);

inbH1 = std(filtfilt(bpChk, cH1(quiet)));
inbL1 = std(filtfilt(bpChk, cL1(quiet)));

% Expected value for correctly whitened data restricted to a band: the
% square root of the fraction of the spectrum that survives the filter.
expBand = sqrt((400-30)/(fsD/2));
fprintf('Std in 30-400 Hz:  H1 = %.3f   L1 = %.3f   (expect ~%.3f)\n', ...
        inbH1, inbL1, expBand);
fprintf('Ratio L1/H1:  full band %.2f,  in band %.2f\n', ...
        stdL1/stdH1, inbL1/inbH1);
fprintf(['(If the in-band ratio is near 1 while the full-band ratio is ' ...
         'not,\n the excess is outside the band the figures show.)\n\n']);

% Renormalize each detector to unit variance in the quiet region before
% the spectrograms. The conditioning already applies sqrt(2/fs), but if
% one detector comes out at a different level the combined map would be
% dominated by it. This is a display choice, applied after conditioning,
% not a change to the chain.
cH1 = cH1 / stdH1;
cL1 = cL1 / stdL1;
fprintf('Renormalized each detector to unit variance for display.\n\n');

% Combined transfer function, as LALdatacondchain computes it 
% The chain pushes a unit impulse through every filter in the chain and
% stores the magnitude of its spectrum as condTF. That is the end-to-end
% frequency response of the conditioning, and it is what lets a later
% stage undo or account for what conditioning did.
%
% Computed here as a diagnostic: the plot should be flat across the
% whitened band and roll off outside it.
impVec = zeros(size(datH1));
impVec(floor(numel(impVec)/2)) = 1;
if dtrndON, impVec = fftfilt(hpB, impVec);
end
if whtnON,  impVec = fftfilt(whtB_H1, impVec);
end
if dnSmplON
    impVec = resample(impVec, 1, P.dnSmplFac, P.dnSmplFltrOrdr);
end
impVec = sqrt(2/fsD) * impVec;

combTF   = fft(impVec);
kNyq     = floor(numel(impVec)/2) + 1;
combTFMag = abs(combTF(1:kNyq));
fTF      = (0:kNyq-1)' * (fsD/numel(impVec));

figure('Name','Conditioning transfer function','Position',[100 100 900 420]);
semilogx(fTF, 20*log10(combTFMag/max(combTFMag)), 'k');
grid on;
xlim([10 fsD/2]); 
ylim([-80 5]);
xlabel('Frequency (Hz)'); 
ylabel('Response (dB, normalized)');
title('End-to-end conditioning transfer function (condTF)');
exportgraphics(gcf, fullfile(figDir,'conditioning_tf.png'), 'Resolution',150);

% Cut to the display window 
% Everything downstream works on the last few seconds. Transforming 700 s
% would be slow and would hide the signal in a very wide plot.
selD  = valid & t >= dsp.tStart-2 & t <= dsp.tStop+2;
% Check the window is FULLY covered, not merely non-empty. A partially
% covered window is what produced figures blank on the right: the data
% simply stopped before the merger.
assert(any(selD), 'Display window falls entirely outside the valid region.');
tValidEnd = t(find(valid,1,'last'));
assert(tValidEnd >= dsp.tStop, ...
       ['Valid data ends at %+.1f s but the display window runs to ' ...
        '%+.1f s. Increase anaStop by at least %.0f s.'], ...
       tValidEnd, dsp.tStop, ceil(dsp.tStop - tValidEnd + 2));
xH1   = cH1(selD);
xL1   = cL1(selD);
tDsp0 = t(find(selD,1));

fprintf('Display segment: %.1f s from t = %+.1f s\n\n', numel(xH1)/fsD, tDsp0);

%  Q-TRANSFORM  
% A Q-transform holds Q = f/df constant, so its window shortens as
% frequency rises. That is what a chirp needs: GW170817 barely moves in
% frequency for most of its length, then sweeps hundreds of Hz in the
% final second. An STFT, with one fixed window, cannot serve both ends.
%
% This is also the transform behind GravitySpy's glitch images.

fprintf('Computing Q-transforms (Q = %g)...\n', qt.Q);
tic;
[QH1, fQ, tQ] = qScan(xH1, fsD, qt.Q, [dsp.fMin dsp.fMax], qt.nFreq);
[QL1, ~,  ~ ] = qScan(xL1, fsD, qt.Q, [dsp.fMin dsp.fMax], qt.nFreq);
tQ_H1 = tQ + tDsp0;
tQ_L1 = tQ + tDsp0 + dtLH;       % L1 arrived 3 ms earlier

% Combine. "Combined" is not defined in the Fig. 2 caption, so this is a
% choice: qScan already normalizes each row to its own background, so
% summing weights the two detectors equally.
QL1onH1 = interp1(tQ_L1, QL1', tQ_H1, 'linear', 0)';
QComb = QH1 + QL1onH1;
fprintf('        done (%.1f s)\n\n', toc);
figure('Name','GW170817 Q-transform','Position',[100 100 1000 800]);

subplot(3,1,1);
imagesc(tQ_H1, fQ, QH1);
axis xy; 
set(gca,'YScale','log');
ylim([dsp.fMin dsp.fMax]); 
xlim([dsp.tStart dsp.tStop]);
clim(climSingle);
colormap(parula);
colorbar;
ylabel('Frequency (Hz)');
title(sprintf('H1 — Hanford   (Q = %g)', qt.Q));

subplot(3,1,2);
imagesc(tQ_L1, fQ, QL1); 
axis xy;
set(gca,'YScale','log');
ylim([dsp.fMin dsp.fMax]);
xlim([dsp.tStart dsp.tStop]);
clim(climSingle);
colormap(parula); 
colorbar;
ylabel('Frequency (Hz)');
title('L1 — Livingston');

subplot(3,1,3);
imagesc(tQ_H1, fQ, QComb);
axis xy;
set(gca,'YScale','log');
ylim([dsp.fMin dsp.fMax]);
xlim([dsp.tStart dsp.tStop]);
clim(climComb);
colormap(parula);
colorbar;
xlabel('Time from merger (s)');
ylabel('Frequency (Hz)');
title('Combined H1 + L1 — compare with the GW inset of Fig. 2');
exportgraphics(gcf, fullfile(figDir,'fig2_qtransform.png'), 'Resolution',150);
%  STFT, for comparison
% Same data, fixed window. Worth producing side by side: the STFT track
% is blurred at the high-frequency end where the Q-transform is sharp.
winLen = round(st.winSec*fsD);
nOver  = round(0.95*winLen);
nFFT   = 2^nextpow2(4*winLen);

[Sh, f, tSh] = spectrogram(xH1, hann(winLen), nOver, nFFT, fsD);
[Sl, ~, tSl] = spectrogram(xL1, hann(winLen), nOver, nFFT, fsD);

tSh = tSh + tDsp0;
tSl = tSl + tDsp0 + dtLH;

keep = f >= dsp.fMin & f <= dsp.fMax;
Ah = abs(Sh(keep,:)); Ah = Ah ./ median(Ah(:));
Al = abs(Sl(keep,:)); Al = Al ./ median(Al(:));
AlOnH = interp1(tSl, Al', tSh, 'linear', 0)';
Acomb = Ah + AlOnH;

figure('Name','GW170817 STFT','Position',[100 100 1000 800]);

subplot(3,1,1);
imagesc(tSh, f(keep), Ah); 
axis xy;
set(gca,'YScale','log');
ylim([dsp.fMin dsp.fMax]); 
xlim([dsp.tStart dsp.tStop]);
clim([0 2.5]); colormap(parula); 
colorbar;
ylabel('Frequency (Hz)');
title('H1 — Hanford (STFT)');

subplot(3,1,2);
imagesc(tSl, f(keep), Al);
axis xy;
set(gca,'YScale','log');
ylim([dsp.fMin dsp.fMax]);
xlim([dsp.tStart dsp.tStop]);
clim([0 2.5]); 
colormap(parula); 
colorbar;
ylabel('Frequency (Hz)');
title('L1 — Livingston (STFT)');

subplot(3,1,3);
imagesc(tSh, f(keep), Acomb);
axis xy;
set(gca,'YScale','log');
ylim([dsp.fMin dsp.fMax]);
xlim([dsp.tStart dsp.tStop]);
clim([0 4]); colormap(parula);
colorbar;
xlabel('Time from merger (s)');
ylabel('Frequency (Hz)');
title('Combined H1 + L1 (STFT, for comparison)');
exportgraphics(gcf, fullfile(figDir,'fig2_stft.png'), 'Resolution',150);
%  Helper functions
function [strain, tRel] = loadStrain(fname, tEvent, fs)
% Reads strain and builds the merger-relative time axis.
% NOTE: Xstart and Xspacing are attributes of /strain/Strain, not /meta.
% Both groups exist in these files; reading the wrong one gives a
% silently incorrect time axis.
    assert(exist(fname,'file')==2, 'File not found: %s', fname);
    strain = h5read(fname,'/strain/Strain');
    strain = strain(:);
    dt     = double(h5readatt(fname,'/strain/Strain','Xspacing'));
    gps0   = double(h5readatt(fname,'/strain/Strain','Xstart'));
    assert(abs(dt - 1/fs) < 1e-12, 'Sample rate mismatch in %s', fname);
    tRel   = gps0 + (0:numel(strain)-1)'*dt - tEvent;
end


function validateDataQuality(fname, strain, label, tRel)
% GWOSC data quality and hardware injection masks.
% NOTE THE INVERTED CONVENTION on the injection mask: the bit names read
% NO_CW_HW_INJ, NO_BURST_HW_INJ and so on, so a CLEARED bit means that
% injection type IS present.
% When tRel is supplied, any second where DQmask < 127 is located on the
% merger-relative axis. That matters because a bad stretch inside the PSD
% training window gives a noise estimate that does not describe the data,
% which shows up later as a conditioned standard deviation far from 1.
    fprintf('%s data quality:\n', label);
    try
        dq  = h5read(fname,'/quality/simple/DQmask');
        inj = h5read(fname,'/quality/injections/Injmask');
        nm  = h5read(fname,'/quality/injections/InjShortnames');
        fprintf('    DQmask  min = %d\n', min(dq));
        fprintf('    Injmask min = %d\n', min(inj));

        % Locate the bad seconds. DQmask carries one flag per second.
        if nargin >= 4 && min(dq) < 127
            tDQ = tRel(1) + (0:numel(dq)-1)';
            bad = tDQ(dq < 127);
            fprintf('    >> %d s with DQmask < 127, spanning t = %+.0f to %+.0f s\n', ...
                    numel(bad), min(bad), max(bad));
            % Report contiguous blocks, so a single dropout is not
            % mistaken for a span covering the whole record.
            d = [1; find(diff(bad) > 1.5)+1];
            e = [find(diff(bad) > 1.5); numel(bad)];
            for j = 1:min(numel(d), 10)
                fprintf('       block %d: t = %+.0f to %+.0f s\n', ...
                        j, bad(d(j)), bad(e(j)));
            end
            if numel(d) > 10
                fprintf('       ... and %d more blocks\n', numel(d)-10);
            end
        end
        mv = double(min(inj));
        any_ = false;
        for b = 1:numel(nm)
            if bitget(mv,b) == 0
                fprintf('    >> active injection: %s\n', strtrim(nm{b}));
                any_ = true;
            end
        end
        if ~any_, fprintf('    no hardware injections active\n');
        end
    catch
        fprintf('    quality channels not present in this file\n');
    end
    n = sum(isnan(strain));
    if n > 0
        warning('%s: %d NaN samples.', label, n);
    else
        fprintf('    no NaN samples\n');
    end
    fprintf('\n');
end

function [pxx, f] = pwelchMedian(x, nSeg, nOvrlp, win, fs)
% Thin wrapper around the repository's pwelch_median.

    if exist('pwelch_median', 'file') == 2
        [pxx, f] = pwelch_median(x, nSeg, nOvrlp, win, fs);
        return;
    end
    %  Fallback
    x    = x(:)';
    win  = win(:)';
    hop  = nSeg - nOvrlp;
    nWin = floor((length(x) - nOvrlp) / hop);
    kNyq = floor(nSeg/2) + 1;
    winNorm = sum(win.^2);
    P = zeros(nWin, kNyq);

    for k = 1:nWin
        iStart = (k-1)*hop + 1;
        iEnd   = iStart + nSeg - 1;
        if iEnd > length(x);
            break; 
        end
        Xk = fft(x(iStart:iEnd) .* win, nSeg);
        Pk = (abs(Xk(1:kNyq)).^2) / (fs * winNorm);
        Pk(2:end-1) = 2 * Pk(2:end-1);
        P(k,:) = Pk;
    end

    pxx = median(P(1:k,:), 1)';
    f   = (0:(kNyq-1))' * (fs/nSeg);
end

function b = designWhitener(pxx, pxxF, fs, P, whtFltrOrdr)
% Designs the whitening FIR filter as LALdatacondchain does.
%
% The key step is LEVELLING. Outside [hiPassFltrCut, whtnHiFreq] the PSD
% is held at its band-edge value rather than used as measured, giving a
% flat finite response there instead of a sharp cut.

    idxLevelLoStrt = find(pxxF <= P.hiPassFltrCut, 1, 'last');
    idxLevelLo     = pxxF <  P.hiPassFltrCut;
    idxLevelHiStrt = find(pxxF >= P.whtnHiFreq, 1);
    idxLevelHi     = pxxF >  P.whtnHiFreq;

    pxxLvld = pxx;
    if ~isempty(idxLevelLoStrt)
        pxxLvld(idxLevelLo) = pxxLvld(idxLevelLoStrt);
    end
    if ~isempty(idxLevelHiStrt)
        pxxLvld(idxLevelHi) = pxxLvld(idxLevelHiStrt);
    end

    whtTF = 1 ./ sqrt(pxxLvld);

    % fir2 needs a normalized frequency vector running from 0 to 1.
    fN = pxxF(:) / (fs/2);
    fN(1) = 0;  fN(end) = 1;

    b = fir2(whtFltrOrdr, fN, whtTF(:));
end
function [Qmap, fQ, tQ] = qScan(x, fs, Q, fRange, nFreq)
% Q-transform of a conditioned time series.
    x = x(:);
    n = numel(x);
    X = fft(x);

    f  = (0:n-1)' * (fs/n);
    fQ = logspace(log10(fRange(1)), log10(fRange(2)), nFreq)';

    Qmap = zeros(nFreq, n);

    for k = 1:nFreq
        f0 = fQ(k);
        hw = (f0/Q) * sqrt(11)/2;

        w   = zeros(n,1);
        sel = abs(f - f0) < hw & f < fs/2;
        xi  = (f(sel) - f0) / hw;
        w(sel) = (1 - xi.^2).^2;

        nrm = sqrt(sum(w.^2));
        if nrm == 0,
            continue;
        end
        w = w / nrm;

        row = ifft(X .* w) * sqrt(n);
        Qmap(k,:) = abs(row).^2;
    end

    med = median(Qmap, 2);
    med(med == 0) = 1;
    Qmap = Qmap ./ med;

    tQ = (0:n-1)' / fs;
end
