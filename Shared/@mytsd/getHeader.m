function csc_info = gh2(tsa)

csc_info = [];
for hline = 1:length(tsa.units)
   
    line = strtrim(tsa.units{hline});
    
    if isempty(line) | ~strcmp(line(1),'-') % not an informative line, skip
        continue;
    end
    
    a = regexp(line(2:end),'(?<key>\w+)\s+(?<val>\S+)','names');
    
    % deal with characters not allowed by MATLAB struct
    if strcmp(a.key,'DspFilterDelay_µs')
        a.key = 'DspFilterDelay_us';
    end
    
    csc_info = setfield(csc_info,a.key,a.val);
    
    % convert to double if possible
    if ~isnan(str2double(a.val))
        csc_info = setfield(csc_info,a.key,str2double(a.val));
    end
    
end