function PhasePowerCorr(fname, varargin)
% function PhasePowerCorr(fname, varargin)
%
% INPUTS
%
% fname: name of file to extract LFPs from. **Note: if you want
% to use sample file name, remove comment marker from fc line**
%
% varargins (with defaults):
% F1 = frequency band to analyze power. Default is [45 65]
% F1_stop = indicates the stopband for filter for F1. Default is [40 70]
% F2 = frequency band to analyze phase. Default is [3 4]
% F2_stop = indicates the stopband for filter for F2. Default is [1 6]
% restrict_data = identify selection of data to use. Default is [2700 3300]
% dsf = down sampling factor. Default is 4
%
% OUTPUTS
%
% A scatterplot depicting the phase power relationship between F1 and F2.
% Returns the mean power for each phase bin.


F1 = [45 65]; % gamma
F1_stop = [40 70];
F2 = [3 4]; % delta
F2_stop = [1 6];
restrict_data = [2700 3300]; %risk session
dsf = 4;

extract_varargin;

%% load and restrict the data
%  fname = 'R016-2012-10-03-CSC04a.Ncs';
csc = LoadCSCnew(fname);
csc_Header = getHeader(csc);
Fs = getfield(csc_Header,'SamplingFrequency');

cscR = Restrict(csc,restrict_data(1),restrict_data(2));

%% decimate data
Fs = Fs/dsf;
d = decimate(Data(cscR),dsf);

%% design filters for frequency ranges
Wp = F1 * 2 / Fs; % gamma
Ws = F1_stop * 2 / Fs;
[N,Wn] = cheb1ord(Wp, Ws, 3, 20);
[b,a] = cheby1(N,0.5,Wn);

Wp = F2 * 2 / Fs; % delta
Ws = F2_stop * 2 / Fs;
[N,Wn] = cheb1ord(Wp, Ws, 3, 20);
[b2,a2] = cheby1(N,0.5,Wn);

%% filter the data
gamma_filtered = filtfilt(b,a,d);
delta_filtered = filtfilt(b2,a2,d);

%% convert to power envelope
delta_power = delta_filtered.^2;
delta_power_filtered = medfilt1(delta_power,25);

%% convert to z-score
delta_zscore = zscore(delta_power_filtered);

%% find times when above threshold
index = delta_zscore > 5;

%% create new variables with only delta-containing epochs
var1 = 1;

for ii = 1:length(index);
    if index(ii) ==1;
        delta_filtered_2(var1) = delta_filtered(ii);
        gamma_filtered_2(var1) = gamma_filtered(ii);
        var1 = var1+1;
    end
end

%% extract delta phase and low gamma power
gamma_amp = abs(hilbert(gamma_filtered)); %original signal
delta_phi = angle(hilbert(delta_filtered));

gamma_amp2 = abs(hilbert(gamma_filtered_2)); %only delta-containing epochs
delta_phi2 = angle(hilbert(delta_filtered_2));

%% use MODIFIED averageXbyYbin to plot relationship
phi_edges = -pi:pi/8:pi;
[pow_bin,pow_bin_sd] = averageXbyYbin(gamma_amp,delta_phi,phi_edges); %original signal
[pow_bin2,pow_bin_sd2] = averageXbyYbin(gamma_amp2,delta_phi2,phi_edges); %only delta-containing epochs

pow_bin(end-1) = pow_bin(end-1)+pow_bin(end);
pow_bin = pow_bin(1:end-1);
pow_bin2(end-1) = pow_bin2(end-1)+pow_bin2(end);
pow_bin2 = pow_bin2(1:end-1);
phi_centers = phi_edges(1:end-1)+pi/16;

% to plot all data within risk session, remove following comment markers
% subplot(121);
% errorbar(phi_centers,pow_bin,pow_bin_sd); xlabel('Phase bins'); ylabel('Power'); title('All within risk session');
% subplot(122);
errorbar(phi_centers,pow_bin2,pow_bin_sd2); xlabel('Phase bins'); ylabel('Power'); title('Phase-power relationship (delta-containing epochs only)');