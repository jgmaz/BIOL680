function neuroplot(S,csc,varargin)
spikecolor = [1 0 0];
LFPcolor = [.7 .7 .7];
evtcolor = [0 0 0];

extract_varargin;

%% get spikes
for st = 1:length(S)
    rspike{st} = Data(S{st});
end

%% Raster plot
figure
hold on
for j = 1:67
    for i = 1:length(rspike{j})
        line([rspike{j}(i) rspike{j}(i)],[j-1 j],'color',spikecolor)
        
    end
end

%% LFP plot
plot(Range(csc),(Data(csc)/200)+75,'color',LFPcolor)
set(gca,'XLim',[5999 6001],'ytick',[])

%% events
line([6000 6000],[0 67],'color',evtcolor)
line([5990 5990],[0 67],'color',evtcolor)
line([6003 6003],[0 67],'color',evtcolor)

end