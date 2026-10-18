code/GW150914.m
%  GW150914 — Fig. 1 time series
%  Reproduces the three time-series rows of Fig. 1 in Abbott et al.,
%  PRL 116, 061102 (2016), and adds a Q-transform of the same data.
dataDir = fullfile(getenv('HOME'),'Desktop','GW150914','data');
figDir  = fullfile(getenv('HOME'),'Desktop','GW150914','figs');
if ~exist(figDir,'dir');
    mkdir(figDir);
end

% The 4096 s files are required: 264 s of training data cannot come from
% the 32 s file.
fileH1 = fullfile(dataDir,'H-H1_LOSC_4_V2-1126257414-4096.hdf5');
fileL1 = fullfile(dataDir,'L-L1_LOSC_4_V2-1126257414-4096.hdf5');
fileNR = fullfile(dataDir,'fig2-unfiltered-template-reconstruction-H.txt');

% Event parameters 
fs    = 4096;
tZero = 1126259462;       % GPS origin of the published figure axis
delay = 0.0069;           % L1 first, H1 6.9 ms later

% Stretch to condition, relative to tZero. Must hold 264 s of training
% plus 18 s of settling, with the signal well clear of the training data.
anaStart = -400;
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

% Display filtering 
D.fLow  = 35;
D.fHigh = 350;
D.notch = [60 120 180 240 300 330];

% Q-transform 
%  Q = 10  the compromise used here
qt.Q      = 10;
qt.fRange = [20 400];
qt.nFreq  = 120;

% Colour ceiling. With the lines notched out the useful dynamic range
% changed, so this was raised from 10. If the map saturates, raise it
% further; if the track is faint, lower it.
qt.clim   = [0 20];

% Run the Q-transform on the NOTCHED data rather than the raw conditioned
% series. The conditioning leaves the 60 Hz harmonics and the 330 Hz
% calibration line in place, and on a normalized-energy map they appear
% as bright horizontal stripes. Set false to see them.
qt.useNotched = true;

% Plot window, matching the published axes
xLo = 0.285;
xHi = 0.46;

%  Load and validate
[rawH1, tRawH1] = loadStrain(fileH1, tZero, fs);
[rawL1, tRawL1] = loadStrain(fileL1, tZero, fs);
assert(abs(tRawH1(1) - tRawL1(1)) < 1e-9, 'H1 and L1 time axes differ');

validateDataQuality(fileH1, rawH1, 'H1', tRawH1);
validateDataQuality(fileL1, rawL1, 'L1', tRawL1);

selAna = tRawH1 >= anaStart & tRawH1 <= anaStop;
datH1  = rawH1(selAna);
datL1  = rawL1(selAna);
tAna0  = tRawH1(find(selAna,1));
nAna   = numel(datH1);
clear rawH1 rawL1 tRawH1 tRawL1

tAna = tAna0 + (0:nAna-1)'/fs;
fprintf('Conditioning %.0f s, t = %+.0f to %+.0f s\n\n', nAna/fs, tAna(1), tAna(end));

% Numerical relativity waveform
NR  = readmatrix(fileNR, 'CommentStyle','#');
tNR = NR(:,1);
hNR = NR(:,2) * 1e-21;

% Taper the first 20 ms. The file covers only 0.25-0.45 s and is zero
% outside, so a hard step would ring through every filter downstream.
nT   = round(0.020*fs);
ramp = 0.5*(1 - cos(pi*(0:nT-1)'/nT));
tp   = ones(size(hNR));  tp(1:nT) = ramp;
hNR  = hNR .* tp;

nrH1 = interp1(tNR, hNR, tAna, 'linear', 0);
% GWOSC publishes the NR waveform projected onto H1 only. This is that
% projection shifted and inverted, not the true L1 antenna projection,
% so L1 agreement is correspondingly worse.
nrL1 = -interp1(tNR, hNR, tAna + delay, 'linear', 0);

%  STEP 1 — DETRENDING
if dtrndON
    % fir1 highpass applied causally with fftfilt. FIR rather than IIR,
    % so the phase is linear and the delay is exactly half the order.
    hpB = fir1(P.hiPassFltrOrdr, P.hiPassFltrCut/(fs/2), 'high');
    tic;
    dtH1   = fftfilt(hpB, datH1);
    dtL1   = fftfilt(hpB, datL1);
    dtNRH1 = fftfilt(hpB, nrH1);
    dtNRL1 = fftfilt(hpB, nrL1);

    fprintf('Step 1  detrending: fir1 highpass %g Hz, order %d  (%.1f s)\n', ...
            P.hiPassFltrCut, P.hiPassFltrOrdr, toc);
    hpDly = P.hiPassFltrOrdr/fs;        % as LALdatacondchain accumulates it
else
    dtH1 = datH1; 
    dtL1 = datL1;  
    dtNRH1 = nrH1; 
    dtNRL1 = nrL1;
    hpDly = 0;
    fprintf('Step 1  detrending: OFF\n');
end

%  STEP 2 — PSD ESTIMATION
% Median-averaged Welch on a training segment. The median rather than the
% mean: each periodogram bin is exponentially distributed about the true
% PSD, so a segment containing a glitch pulls the mean up but 
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

assert(trnLen >= 0.90*needLen, ...
       'Training segment truncated. Move anaStart earlier.');

% The signal must lie outside the training data, or the estimator reads
% the signal's own power as noise and divides it out.
assert(trnSpan(2) < -5, ...
       'Training segment reaches the event. Move anaStart earlier.');
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
grid on; xlim([10 2000]);
ylim([1e-24 1e-19]);
xlabel('Frequency (Hz)');
ylabel('ASD (strain/\surdHz)');
legend('H1','L1');
title('Median-averaged ASD — compare with Fig. 3, inset (b)');
exportgraphics(gcf, fullfile(figDir,'asd.png'), 'Resolution',150);

%  STEP 3 — WHITENING
% The PSD is LEVELLED outside [hiPassFltrCut, whtnHiFreq]: held at the
% band-edge value rather than used as measured. Below 25 Hz the data is
% already gone, so 1/sqrt(PSD) would explode; above 1650 Hz there is
% nothing worth whitening. Levelling gives a flat finite response in both
% regions and keeps the fir2 design well conditioned.
% Data and NR waveform go through the SAME filter. Otherwise the residual
% in row 3 would show the difference between two filters rather than any
% model-data mismatch.
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

    wH1   = fftfilt(whtB_H1, dtH1);
    wL1   = fftfilt(whtB_L1, dtL1);
    wNRH1 = fftfilt(whtB_H1, dtNRH1);
    wNRL1 = fftfilt(whtB_L1, dtNRL1);

    fprintf('      done (%.1f s)\n', toc);
    whtDly = whtFltrOrdr/fs;
else
    wH1 = dtH1;  wL1 = dtL1;  wNRH1 = dtNRH1;  wNRL1 = dtNRL1;
    whtDly = 0;
    fprintf('Step 3  whitening: OFF\n');
end

%  STEP 4 — DOWNSAMPLE
if dnSmplON
    fsD = fs / P.dnSmplFac;

    cH1   = resample(wH1,   1, P.dnSmplFac, P.dnSmplFltrOrdr);
    cL1   = resample(wL1,   1, P.dnSmplFac, P.dnSmplFltrOrdr);
    cNRH1 = resample(wNRH1, 1, P.dnSmplFac, P.dnSmplFltrOrdr);
    cNRL1 = resample(wNRL1, 1, P.dnSmplFac, P.dnSmplFltrOrdr);

    % LALdatacondchain's normalization: makes data at the noise floor
    % have unit variance.
    nrmFac = sqrt(2/fsD);
    cH1   = nrmFac * cH1;
    cL1   = nrmFac * cL1;
    cNRH1 = nrmFac * cNRH1;
    cNRL1 = nrmFac * cNRL1;

    fprintf('Step 4  downsampled %g -> %g Hz (filter order %d), scaled by %.4g\n\n', ...
            fs, fsD, P.dnSmplFltrOrdr, nrmFac);
else
    fsD = fs;
    nrmFac = sqrt(2/fsD);
    cH1 = nrmFac*wH1;  cL1 = nrmFac*wL1;
    cNRH1 = nrmFac*wNRH1;  cNRL1 = nrmFac*wNRL1;
    fprintf('Step 4  downsampling: OFF\n\n');
end

% Time axis, with the group delay removed 
% fftfilt is causal, so both FIR filters delay the output. For a
% linear-phase FIR of order N the group delay is exactly N/2 samples.
% resample() compensates its own delay internally.
%
% LALdatacondchain accumulates N/fs per filter — twice the group delay —
% and uses that as a conservative settling time, recording it as
% fltrTShft rather than correcting the timing.
grpDly  = (P.hiPassFltrOrdr*dtrndON + whtFltrOrdr*whtnON) / (2*fs);
settleT = (P.hiPassFltrOrdr*dtrndON + whtFltrOrdr*whtnON) / fs;

% Data windowing, as LALdatacondchain does
% The chain zeroes everything before the accumulated filter delay rather
% than trimming it, and records that delay as fltrTShft in the output
% file so a downstream consumer knows the data is shifted:
%
%     dtrndWhtdDnsmpl(1:(strtTCutIndx-1)) = 0;
%
% Reproduced here. It removes exactly the stretch the chain considers
% contaminated by filter settling.
strtTCutIdx = min(numel(cH1), floor(max(0,settleT)*fsD) + 1);
nZeroed     = strtTCutIdx - 1;
cH1(1:nZeroed) = 0;
cL1(1:nZeroed) = 0;
cNRH1(1:nZeroed) = 0;
cNRL1(1:nZeroed) = 0;

t = tAna0 + (0:numel(cH1)-1)'/fsD - grpDly;
fprintf('Group delay removed: %.3f s   (settling time %.3f s)\n', grpDly, settleT);

% Valid region: past the settling time at the start, and clear of the
% filter running out at the end.
valid = t > (tAna0 + settleT - grpDly) & t < (t(end) - settleT);
fprintf('Zeroed the first %d samples (%.1f s), as the chain does.\n', ...
        nZeroed, nZeroed/fsD);
fprintf('Valid span: t = %+.1f to %+.1f s\n\n', ...
        t(find(valid,1)), t(find(valid,1,'last')));

assert(any(valid & t > xLo & t < xHi), ...
       'The plot window falls outside the valid region.');

% Conditioning sanity check 
% Conditioned data should have standard deviation near 1 away from the
% signal. Far from 1 means the PSD does not describe this stretch.
quiet = valid & t < -5;
fprintf('Conditioned std (quiet region):  H1 = %.3f   L1 = %.3f   (expect ~1)\n\n', ...
        std(cH1(quiet)), std(cL1(quiet)));


% Combined transfer function, as LALdatacondchain computes it
% The chain pushes a unit impulse through every filter in the chain and
% stores the magnitude of its spectrum as condTF. That is the end-to-end
% frequency response of the conditioning, and it is what lets a later
% stage undo or account for what conditioning did.
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

%  FIGURE 1
bp = designfilt('bandpassiir','FilterOrder',8, ...
                'HalfPowerFrequency1',D.fLow, ...
                'HalfPowerFrequency2',D.fHigh, ...
                'SampleRate',fsD);

fH1   = applyNotches(filtfilt(bp, cH1),   D.notch, fsD);
fL1   = applyNotches(filtfilt(bp, cL1),   D.notch, fsD);
fNRH1 = applyNotches(filtfilt(bp, cNRH1), D.notch, fsD);
fNRL1 = applyNotches(filtfilt(bp, cNRL1), D.notch, fsD);

% Conditioned data is in units of the noise standard deviation, not
% strain. One factor scales everything so the NR waveform peaks at 1.0 —
% the peak strain the paper quotes. The calibration comes from the model,
% not from the data. Nothing is fitted.
scale = 1.0 / max(abs(fNRH1(valid)));
fH1   = fH1   * scale;
fL1   = fL1   * scale;
fNRH1 = fNRH1 * scale;
fNRL1 = fNRL1 * scale;

show = valid & t > xLo-0.02 & t < xHi+0.02;

figure('Name','GW150914 Fig. 1 time series','Position',[80 80 1000 820]);

% ROW 1 — strain
subplot(3,2,1);
plot(t(show), fH1(show), 'r', 'LineWidth',1.1);
grid on;
xlim([xLo xHi]);
ylim([-1.5 1.5]);
ylabel('Strain (10^{-21})');
title('Hanford, Washington (H1)');
legend('H1 observed','Location','southwest');

subplot(3,2,2);
plot(t(show), fL1(show), 'b', 'LineWidth',1.1); 
hold on
% Shift on the time axis, not with circshift, which would wrap the array.
plot(t(show)-delay, -fH1(show), 'Color',[1 .5 .3], 'LineWidth',1.1);
grid on; 
xlim([xLo xHi]); 
ylim([-1.5 1.5]);
title('Livingston, Louisiana (L1)');
legend('L1 observed','H1 observed (shifted, inverted)','Location','southwest');

% ROW 2 — NR waveform over the data
% The published row 2 also shows shaded 90% credible bands from two
% independent reconstructions. Those are Bayesian posteriors, not
% available from the released data, so only the NR waveform is drawn.
subplot(3,2,3);
plot(t(show), fH1(show), 'r', 'LineWidth',1.1); 
hold on
plot(t(show), fNRH1(show), 'k', 'LineWidth',1.4);
grid on;
xlim([xLo xHi]);
ylim([-1.5 1.5]);
ylabel('Strain (10^{-21})');
legend('H1 data','Numerical relativity','Location','southwest');

subplot(3,2,4);
plot(t(show), fL1(show), 'b', 'LineWidth',1.1);
hold on
plot(t(show), fNRL1(show), 'k', 'LineWidth',1.4);
grid on;
xlim([xLo xHi]);
ylim([-1.5 1.5]);
legend('L1 data','Numerical relativity','Location','southwest');

% ROW 3 — residuals
resH1 = fH1 - fNRH1;
resL1 = fL1 - fNRL1;

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

% Residual RMS test 
% NOT part of the published figure. Row 3 makes a claim by eye — that
% once the model is removed nothing structured is left. These are the
% numbers behind it.
inSignal = valid & t > 0.35 & t < 0.45;
inNoise  = valid & t > 0.28 & t < 0.33;

%  Q-TRANSFORM
% The published Fig. 1 bottom row uses a Q-transform (Chatterji et al.,
% CQG 21, S1809, 2004), not a short-time Fourier transform. An STFT uses
% one window length everywhere, so time and frequency resolution are
% fixed across the plane. A Q-transform holds Q = f/df constant, so the
% window shortens as frequency rises — which is what a chirp needs, since
% it sweeps fastest at the end.
% This is also the transform behind GravitySpy's glitch images.

qSel = valid & t > 0.12 & t < 0.62;

figure('Name','GW150914 Q-transform','Position',[80 80 1000 700]);

% fH1 and fL1 were scaled to strain units; undo that so the Q-transform
% sees the same amplitude convention either way. The per-row median
% normalization inside qScan makes the overall scale irrelevant, but
% keeping it consistent avoids confusion.
if qt.useNotched
    qSrcH1 = fH1/scale;  qSrcL1 = fL1/scale;  qTag = 'notched';
else
    qSrcH1 = cH1;        qSrcL1 = cL1;        qTag = 'un-notched';
end

for d = 1:2
    if d == 1
        x = qSrcH1(qSel);  
        name = 'Hanford, Washington (H1)';
    else
        x = qSrcL1(qSel); 
        name = 'Livingston, Louisiana (L1)';
    end

    [Qmap, fQ, tQ] = qScan(x, fsD, qt.Q, qt.fRange, qt.nFreq);
    tQ = tQ + t(find(qSel,1));

    subplot(2,1,d);
    imagesc(tQ, fQ, Qmap);
    axis xy;
    set(gca,'YScale','log');
    ylim(qt.fRange);
    xlim([xLo xHi]);
    clim(qt.clim);
    colormap(parula);
    colorbar;
    ylabel('Frequency (Hz)');
    title(sprintf('%s   (Q = %g, %s)', name, qt.Q, qTag));
    if d == 2,
        xlabel('Time (s)');
    end
end
exportgraphics(gcf, fullfile(figDir,'fig1_qtransform.png'), 'Resolution',150);

%  Helper functions
function [strain, tRel] = loadStrain(fname, tZero, fs)
% Reads strain and builds the time axis from file metadata.
% Xstart and Xspacing are attributes of /strain/Strain, not /meta.
% Both groups exist in these files; reading the wrong one gives a
% silently incorrect time axis.
    assert(exist(fname,'file')==2, 'File not found: %s', fname);
    strain = h5read(fname,'/strain/Strain');
    strain = strain(:);
    dt     = double(h5readatt(fname,'/strain/Strain','Xspacing'));
    gps0   = double(h5readatt(fname,'/strain/Strain','Xstart'));
    assert(abs(dt - 1/fs) < 1e-12, 'Sample rate mismatch in %s', fname);
    tRel   = gps0 + (0:numel(strain)-1)'*dt - tZero;
end


function validateDataQuality(fname, strain, label, tRel)
% GWOSC data quality and hardware injection masks.
%
% NOTE THE INVERTED CONVENTION on the injection mask: the bit names read
% NO_CW_HW_INJ, NO_BURST_HW_INJ and so on, so a CLEARED bit means that
% injection type IS present.
%
% When tRel is supplied, any second where DQmask < 127 is located on the
% figure's time axis. That matters because a bad stretch inside the PSD
% training window gives a noise estimate that does not describe the data,
% which shows up later as a conditioned standard deviation far from 1.
    fprintf('%s data quality:\n', label);
    try
        dq  = h5read(fname,'/quality/simple/DQmask');
        inj = h5read(fname,'/quality/injections/Injmask');
        nm  = h5read(fname,'/quality/injections/InjShortnames');
        fprintf('    DQmask  min = %d   (127 = all categories pass)\n', min(dq));
        fprintf('    Injmask min = %d\n', min(inj));
        if min(dq) < 127
            warning('%s: data quality flags are not all set.', label);
        end

        % Locate the bad seconds. DQmask carries one flag per second.
        if nargin >= 4 && min(dq) < 127
            tDQ = tRel(1) + (0:numel(dq)-1)';
            bad = tDQ(dq < 127);
            fprintf('    >> %d s with DQmask < 127, spanning t = %+.0f to %+.0f s\n', ...
                    numel(bad), min(bad), max(bad));
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
        mv = double(min(inj)); any_ = false;
        for b = 1:numel(nm)
            if bitget(mv,b) == 0
                fprintf('    >> active injection: %s\n', strtrim(nm{b}));
                any_ = true;
            end
        end
        if ~any_, fprintf('  no hardware injections active\n'); 
        end
    catch
        fprintf('  quality channels not present in this file\n');
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
%
% If pwelch_median.m is on the MATLAB path it is called directly, so the
% PSD is computed by exactly the same code LALdatacondchain uses. If it
% is not found, the local fallback below reproduces it line for line.

    if exist('pwelch_median', 'file') == 2
        [pxx, f] = pwelch_median(x, nSeg, nOvrlp, win, fs);
        return;
    end

    % Fallback: a line-for-line copy of pwelch_median.m 
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


function y = applyNotches(x, notchHz, fs)
% Band-reject filters on the instrumental lines, zero-phase.
    y = x;
    for k = 1:numel(notchHz)
        f0 = notchHz(k);
        if f0 >= fs/2 - 2, 
            continue;
        end
        nf = designfilt('bandstopiir','FilterOrder',4, ...
                        'HalfPowerFrequency1',f0-1.5, ...
                        'HalfPowerFrequency2',f0+1.5, ...
                        'SampleRate',fs);
        y = filtfilt(nf, y);
    end
end


function [Qmap, fQ, tQ] = qScan(x, fs, Q, fRange, nFreq)
% Q-transform of a conditioned time series.
% Each row is normalized by its own median, so the colour scale reads as
% "how many times the typical background" — the standard normalized-energy
% display, and what GravitySpy images show.

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
