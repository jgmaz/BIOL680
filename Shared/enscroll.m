%% Enable scroll and zoom
function enscroll
set(gcf,'KeyPressFcn',@enscroll);

function enscroll(~,event)
if strcmp(event.Key,'rightarrow');
    limit = get(gca,'XLim');
    half_range = (limit(2)-limit(1))/2;
    set(gca,'XLim',[limit(1)+half_range limit(2)+half_range]);
elseif strcmp(event.Key,'leftarrow');
    limit = get(gca,'XLim');
    half_range = (limit(2)-limit(1))/2;
    set(gca,'XLim',[limit(1)-half_range limit(2)-half_range]);
elseif strcmp(event.Key,'uparrow');
    zoom(2);
elseif strcmp(event.Key,'downarrow');
    zoom(.5);
end
end
end