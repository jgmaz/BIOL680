function eventLFPplot(fc,varargin)
% function eventLFPplot(fc,varargin)
%
% INPUTS
%
% fc: {1 x 1} cell, name of file to extract LFPs from. **Note: if you want
% to use sample file name, remove comment marker from fc line**
%
% varargins (with defaults):
%
% t_window: [2 x 1] double indicating time window to use, e.g. [1 3] for 1
% second before to 3 seconds after event times. Default is [1 3]
%
% eventtype = 'string' indicating which events to look for, could be
% 'nosepoke', 'reward', 'cue'. Default is 'nosepoke'
%
% location = 'string' indicating which direction the rat took in a
% particular trial, could be 'left', 'right', 'both'. Default is 'both'
%
% block = 'string' indicating which trial type to use, could be 'value',
% 'risk', 'both'. Default is 'both'
%
% cue = cell array with choice of elements {'c1','c3','c5','lo','hi'}.
% Default is {'c1','c3','c5'}
%
% filter_data_pass = [2 x 1] indicating whether or not to apply a Chebyshev type
% 1 filter to the data, e.g. [50 80] sets the passband to 50-80 Hz. Default
% is 0 which = no filter
%
% filter_data_stop = [1 x 1] indicates the stopband, relative to passband,
% for a Chebyshev type 1 filter, e.g. 5 sets the stopband for 5 below and 5
% above the passband, so if passband is 50-80 Hz, then stopband is 45-85
% Hz. Default is 5. **Note: filter_data_pass must not be 0 for a Chebyshev
% type 1 filter to be applied**
%
% decimate_data = [1 x 1] indicating whether or not to decimate data,
% data is downsampled by a factor of x (your input) e.g. 4 would downsample
% the data by a factor of 4. Default is 0.
%
% LFPcolor = modifies the color of LFPs. Default is [.7 .7 .7]

% evtcolor = modifies the color of the trial start time event. Default is
% [0 0 0]

%% varargins

t_window = [1 3];
eventtype = 'nosepoke';
location = 'both';
block = 'both';
cue = {'c1','c3','c5'};
filter_data_pass = 0;
filter_data_stop = 5;
decimate_data = 0;
LFPcolor = [.7 .7 .7];
evtcolor = [0 0 0];
extract_varargin;

%% extract the corresponding piece of LFP

% fc = {'R016-2012-10-03-CSC04a.ncs'}; %file name used for example
data = ft_read_neuralynx_interp(fc);

%% devide the data by trials

cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = t_window(1);
cfg.trialdef.post = t_window(2);

cfg.trialdef.eventtype = eventtype;
cfg.trialdef.location = location;
cfg.trialdef.block = block;
cfg.trialdef.cue = cue;

[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;

data_trl = ft_redefinetrial(cfg,data);

%% optional filter

if filter_data_pass ~= 0;
    Wp = filter_data_pass * 2 / data_trl.hdr.Fs;
    Ws = [filter_data_pass(1)-filter_data_stop filter_data_pass(2)+filter_data_stop] * 2 / data_trl.hdr.Fs;
    [N,Wn] = cheb1ord(Wp, Ws, 3, 20);
    [b,a] = cheby1(N,0.5,Wn);
    
    for ii = 1:length(data_trl.trial);
        data_trl.trial{ii} = filtfilt(b,a,data_trl.trial{ii});
    end
end

%% optional decimate

if decimate_data > 0;
    for ii = 1:length(data_trl.trial);
        data_trl.trial{ii} = decimate(data_trl.trial{ii},decimate_data);
        data_trl.time{ii} = decimate(data_trl.time{ii},decimate_data);
    end
end

%% plot the LFPs trial by trial

for ii = 1:length(data_trl.trial);
    plot(data_trl.time{1},data_trl.trial{ii}+(1000*ii), 'color', LFPcolor);
    hold on;
end

set(gca,'ytick',[]);
line([0 0],[0 1005*ii],'color',evtcolor);