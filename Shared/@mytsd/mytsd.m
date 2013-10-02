function tsa = mytsd(t, qData,tUnits,Header)

switch nargin
    
    case 0
        tsa.t = [];
        tsa.data = [];
         tsa.header = [];
        tsa.units= [];
       
        
        
    case 1
        if isa(t, 'mytsd')
            tsa = t;
            return;
        elseif isa(t, 'ctsd')
            tsa.t = Range(t);
            tsa.data = Data(t);
            if ~isfield(struct(t), 'units')
                warning('MYTSD:Units','units not specified' )
            else
                tsa.units=Units(t);
            end
        elseif isa(t,'struct')
            tsa.t=t.t;
            tsa.data=t.data;
            if ~isfield(t, 'units')
                warning('MYTSD:Units','units not specified' )
            else
                tsa.units=t.units;
            end
            
        else
            error('Unknown copy-from object');
        end
        
    case 2
        tsa.t = t;
        tsa.data = qData;
        tsa.units= 'sec';
        
    case 3
        tsa.t = t;
        tsa.data = qData;
        tsa.units= tUnits;
        
    case 4
        tsa.t = t;
        tsa.data = qData;
        tsa.units = tUnits;
        tsa.header= Header;
        
    otherwise
        error('Constructor error mytsd');
end

tsa = class(tsa, 'mytsd');

