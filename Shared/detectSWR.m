function evt = detectSWR(csc,varargin)
%% extract varargins
ripple_band = [140 180]; %passband
threshold = 5; %threshold to count as SWR event
range_stop = 2.5; %how far away stopband is from passband (Hz) 
fname = 'R042-2013-08-18-CSC03a.ncs'; %name of file to be loaded
time_for_restrict = 900; %length of window to analyze (sec)

extract_varargin

%% load signal
csc = LoadCSCnew(fname);
csc_Header = getHeader(csc);
Fs = getfield(csc_Header,'SamplingFrequency');
csc_length = Range(csc);
csc_stop = csc_length(end);
csc_start = csc_stop - time_for_restrict;
csc_restrict = Restrict(csc,csc_start,csc_stop); %use last 15 mins

%% filter in the ripple band
Wp = ripple_band * 2 / Fs;
Ws = [ripple_band(1)-range_stop ripple_band(2)+range_stop] * 2 / Fs;
[N,Wn] = cheb1ord(Wp, Ws, 3, 20);
[b,a] = cheby1(N,0.5,Wn);

% display filter visualization tool
prompt = 'Look at filter via fvtool? Y/N:   ';
str = input(prompt,'s');
if str == 'y';
    fvtool(b,a);
    prompt = 'Is filter sufficient to continue? Y/N:   ';
    stri = input(prompt,'s');
    if stri == 'n';
        return
    end
end

csc_Filtered = filtfilt(b,a,Data(csc_restrict));

%% convert to power envelope
csc_power = csc_Filtered.^2;
csc_power_filtered = medfilt1(csc_power,25);

%% convert to z-score to deal with variable baseline
csc_zscore = zscore(csc_power_filtered);

%% find times when above threshold
index = csc_zscore > threshold;

%% find crossings from below to above threshold and vice versa
crossings = diff(index);
crossings = [0; crossings];

%% get center time and power
SWR_start = crossings > 0;
SWR_end = crossings < 0;
var1 = 1;
var2 = 1;

for ii = 1:length(SWR_start)
    if SWR_start(ii) ==1
        SWR2_start(var1) = ii;
        var1 = var1+1;
    end
    if SWR_end(ii) ==1
        SWR2_end(var2) = ii;
        var2 = var2+1;
    end
end


for ii = 1:length(SWR2_start)
    SWR_times(ii) = round(median(SWR2_start(ii):SWR2_end(ii)));
end

for ii = 1:length(SWR_times)
    SWR_power(ii) = csc_zscore(SWR_times(ii));
end

SWR_times = SWR_times * 1/ Fs + csc_start;

%% create evt struct
evt.t = SWR_times;
evt.pwr = SWR_power;

%% load spikes

%prompt to plot data
prompt = 'Plot data? Y/N:   ';
stri = input(prompt,'s');
if stri == 'n';
    return
end

Spikes = LoadSpikes(FindFiles('*.t'));
for ii = 1:length(Spikes);
    S{ii} = Restrict(Spikes{ii},csc_start,csc_stop);
end

%% Neuroplot
spikecolor = [1 0 0];
LFPcolor = [.7 .7 .7];
evtcolor = [0 1 0];
SWR_start_plot = SWR2_start * 1/ Fs + csc_start;
SWR_end_plot = SWR2_end * 1/ Fs + csc_start;

%% get spikes
for st = 1:length(S);
    rspike{st} = Data(S{st});
end

%% Raster plot
figure
hold on
for j = 1:length(S);
    for i = 1:length(rspike{j});
        line([rspike{j}(i) rspike{j}(i)],[j-1 j],'color',spikecolor);
        
    end
end

%% LFP plot
plot(Range(csc_restrict),Data(csc_restrict)/500000+length(S)+10,'color',LFPcolor);
set(gca,'XLim',[5987.5 5992.5],'ytick',[]);

%% events
for ii = 1:length(SWR_times)
    
    line([SWR_start_plot(ii) SWR_start_plot(ii)],[0 67],'color',evtcolor);
    line([SWR_end_plot(ii) SWR_end_plot(ii)],[0 67],'color',evtcolor);
end

%% enable enscroll
enscroll;