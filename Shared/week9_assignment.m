%% for each of the two neurons, restrict the data to [3200 5650] 
% (the time interval when the rat was running on the track)

cd('D:\Jimmie\GitHub\BIOL680\Data\Promoted\R042\R042-2013-08-18')

fc = FindFiles('*.t');
S = LoadSpikes(fc);

cell1_id = 5; cell2_id = 42;
 
s1 = Data(Restrict(S{cell1_id},3200,5650));
s2 = Data(Restrict(S{cell2_id},3200,5650));

%% compute the spike density function for each, making sure that your tvec 
% runs from 3200 to 5650 also, and that you have a 50ms SD for the Gaussian convolution kernel
 
binsize = 0.001; % in seconds, so everything else should be seconds too
t = [3200 5650];

tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
 
s1_count = histc(s1,tbin_edges);
s2_count = histc(s2,tbin_edges);

s1_count = s1_count(1:end-1);
s2_count = s2_count(1:end-1);

gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.05./binsize; % 0.05 seconds (50ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf_s1 = conv2(s1_count,gk,'same'); % convolve with gaussian window
gau_sdf_s2 = conv2(s2_count,gk,'same'); % convolve with gaussian window

figure;
subplot(211);
plot(tbin_centers,gau_sdf_s1);
subplot(212);
plot(tbin_centers,gau_sdf_s2);

%% to use these SDFs to generate Poisson spike trains, convert the firing rates 
% given by the SDF to a probability of emitting a spike in a given bin. 
%(As you did above for a 0.47 Hz constant firing rate.)
% generate Poisson spike trains, making sure to use the same tvec

pspike_s1 = gau_sdf_s1.*10^-3;
pspike_s2 = gau_sdf_s2.*10^-3;

pspike_s1 = pspike_s1';
pspike_s2 = pspike_s2';

dt = 0.001;
tvec = t(1)+.001:dt:t(2);

rng default; % reset random number generator to reproducible state
spk_poiss_s1 = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx_s1 = find(spk_poiss_s1 < pspike_s1); % index of bins with spike
spk_poiss_t_s1 = tvec(spk_poiss_idx_s1)'; % use idxs to get corresponding spike time

rng default; % reset random number generator to reproducible state
spk_poiss_s2 = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx_s2 = find(spk_poiss_s2 < pspike_s2); % index of bins with spike
spk_poiss_t_s2 = tvec(spk_poiss_idx_s2)'; % use idxs to get corresponding spike time


 
%% convert Poisson spike trains to ts objects and compute the ccf
ts_s1 = ts(spk_poiss_t_s1);
ts_s2 = ts(spk_poiss_t_s2);
[xcorr,xbin] = ccf(ts_s1,ts_s2,0.01,1); % for poisson spike trains

ts_s1_org = ts(s1);
ts_s2_org = ts(s2);
[xcorr_original,xbin_original] = ccf(ts_s1_org,ts_s2_org,0.01,1); % for original data

figure;
subplot(211);
plot(xbin_original,xcorr_original); title('Original');
subplot(212);
plot(xbin,xcorr); title('Poisson spike trains');