%% 1. white noise PSD
%it is not 1/f. It doesn't vary across f, it is static.

white_noise = rand(168960,1);
[Pxx2,F2] = periodogram(white_noise,hamming(length(white_noise)),length(white_noise),300);
plot(F2,10*log10(Pxx2),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title ('White Noise');
xlim([0 150]);

%% 2. LFPs PSD 
%in the hippocampal LFP, theta (~7 Hz) is pronounced. 
%That's the only oscillation I see, except for a possible blip at ~22 Hz.

%Load vStr data
csc = LoadCSC('R016-2012-10-08-CSC04d.ncs');
run(FindFile('*keys.m'));

%restrict csc
csc_pre = Restrict(csc,0,ExpKeys.TimeOnTrack(1)-10);

%extract values
csc_preR = Range(csc_pre);
csc_preD = Data(csc_pre);
Fs = 1./mean(diff(csc_preR));

%downsample
dsf = 4;
csc_preD = decimate(csc_preD,dsf);
Fs = Fs./dsf;

%create vStr PSD
[Pxx,F] = periodogram(csc_preD,hamming(length(csc_preD)),length(csc_preD),Fs);
figure
subplot(211);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title ('Ventral Striatum');
xlim([0 150]);

% Load hippocampal data
hpc_csc = LoadCSC('R016-2012-10-08-CSC02b.ncs');

% Restrict csc
hpc_csc_pre = Restrict(hpc_csc,0,ExpKeys.TimeOnTrack(1)-10);

%extract values
hpc_csc_preR = Range(hpc_csc_pre);
hpc_csc_preD = Data(hpc_csc_pre);
Fs_hpc = 1./mean(diff(hpc_csc_preR));

% downsample
hpc_csc_preD = decimate(hpc_csc_preD,dsf);
Fs_hpc = Fs_hpc./dsf;

% create hippocampal PSD
[Pxx_hpc,F_hpc] = periodogram(hpc_csc_preD,hamming(length(hpc_csc_preD)),length(hpc_csc_preD),Fs_hpc);
subplot(212);
plot(F_hpc,10*log10(Pxx_hpc),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title('Hippocampus');
xlim([0 150]);

%% 3. messing around with window size parameter for the Welch power spectrum
%decreasing the window size smoothes out the PSD, making it drastically
%easier to see the delta, theta, beta, and gamma peaks in the vStr PSD, and
%the theta peak in the hippocampal PSD.

wSize = 50000;
figure

%create PSDs
subplot(421); %original vStr PSD
plot(F,10*log10(Pxx),'k'); ylabel('Power (dB)'); title('Ventral Striatum');
xlim([0 150]);

subplot(423); %Welch wSize=50000 vStr PSD
[Pxx3,F3] = pwelch(csc_preD,wSize,wSize/2,length(csc_preD),Fs);
plot(F3,10*log10(Pxx3),'g'); ylabel('Power (dB)'); 
xlim([0 150]);

subplot(425); %Welch wSize=5000 vStr PSD
[Pxx4,F4] = pwelch(csc_preD,wSize/10,wSize/20,length(csc_preD),Fs);
plot(F4,10*log10(Pxx4)); ylabel('Power (dB)'); 
xlim([0 150]);

subplot(427); %Welch wSize=500 vStr PSD
[Pxx5,F5] = pwelch(csc_preD,wSize/100,wSize/200,length(csc_preD),Fs);
plot(F5,10*log10(Pxx5),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 
xlim([0 150]);

subplot(422); %original hippocampal PSD
plot(F_hpc,10*log10(Pxx_hpc),'k'); title('Hippocampus');
xlim([0 150]);

subplot(424); %Welch wSize=50000 hippocampal PSD
[Pxx_hpc2,F_hpc2] = pwelch(hpc_csc_preD,wSize,wSize/2,length(hpc_csc_preD),Fs_hpc);
plot(F_hpc2,10*log10(Pxx_hpc2),'g'); 
xlim([0 150]);

subplot(426); %Welch wSize=5000 hippocampal PSD
[Pxx_hpc3,F_hpc3] = pwelch(hpc_csc_preD,wSize/10,wSize/20,length(hpc_csc_preD),Fs_hpc);
plot(F_hpc3,10*log10(Pxx_hpc3));
xlim([0 150])

subplot(428); %Welch wSize=500 hippocampal PSD
[Pxx_hpc4,F_hpc4] = pwelch(hpc_csc_preD,wSize/100,wSize/200,length(hpc_csc_preD),Fs_hpc);
plot(F_hpc4,10*log10(Pxx_hpc4),'r'); xlabel('Frequency (Hz)');
xlim([0 150]);

%% 4. Downsampling vs. decimate
%downsampling versus using decimate distorts the PSD. The original peaks
%are gone, and a new ~30 Hz peak is introduced. Is that the original theta
%signal?

%downsample
csc_preD2 = downsample(csc_preD,dsf);
hpc_csc_preD2 = downsample(hpc_csc_preD,dsf);

%create PSDs
figure
subplot(221); %original vStr PSD
plot(F,10*log10(Pxx),'k'); ylabel('Power (dB)'); title ('Ventral Striatum');
xlim([0 150]);

[Pxx6,F6] = periodogram(csc_preD2,hamming(length(csc_preD2)),length(csc_preD2),Fs);
subplot(223); %downsampled vStr PSD
plot(F6,10*log10(Pxx6)); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);

subplot(222); %original hippocampal PSD
plot(F_hpc,10*log10(Pxx_hpc),'k');  title ('Hippocampus');
xlim([0 150]);

[Pxx_hpc5,F_hpc5] = periodogram(hpc_csc_preD2,hamming(length(hpc_csc_preD2)),length(hpc_csc_preD2),Fs_hpc);
subplot(224); %downsampled hippocampal PSD
plot(F_hpc5,10*log10(Pxx_hpc5)); xlabel('Frequency (Hz)'); 
xlim([0 150]);