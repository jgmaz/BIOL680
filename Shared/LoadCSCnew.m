function csc = LoadCSCnew(fname,varargin)

SubSmp = 2000; %subsampling rate
ADtoVolts = .45776;  %converts AD units to microvolts
VoltCNV = 10^3; %converts microvolts to millivolts
TimeCNV = 10^-6; %converts nerualynx units to seconds

extract_varargin
%% load data
[Timestamps, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 0 1 1 1], 1, 1, []);

%% unwrap samples
csc_data = reshape(Samples,[size(Samples,1)*size(Samples,2) 1]);
csc_data = csc_data.*ADtoVolts.*VoltCNV;

%% creating matching timestamps variable

csc_timestamps = repmat(Timestamps,[size(Samples,1) 1]).*TimeCNV;

dtvec = (0:size(Samples,1)-1)*(1/SubSmp);
dtmat = repmat(dtvec',[1 size(Samples,2)]);

csc_timestamps = csc_timestamps+dtmat;
csc_timestamps = reshape(csc_timestamps,[size(csc_timestamps,1)*size(csc_timestamps,2) 1]);

%% exclude invalid samples
index = (NumberOfValidSamples == 512);
index = repmat(index,[size(Samples,1) 1]);
index = reshape(index,[size(index,1)*size(index,2) 1]);
index = ~(index-csc_data);
csc_data(index) = [];
csc_timestamps(index) = [];

%% create tsd
csc = mytsd(csc_timestamps,csc_data,Header);