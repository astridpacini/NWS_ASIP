function [code,conc,rang] = egg_code_2007_to_2015(ct,form)
%PARSE egg code here. This uses the onfomration in NSIDC usergguide to
%assign concentration intervals to an egg code value, and then convert that
%to a sea ice concentration (in tenths). Note that this function handles
%the pre Oct 1 2015 NWS ASIP datafiles, which tabulate concentration
%information with a character string with a dash (e.g. '9-10'). 
%Written by A. Pacini, apacini@uw.edu, May 2024

if strcmp(ct,'') && strcmp(form,'OPEN') % open water
    code = 0;
    conc = 0;
    rang = 0;
elseif strcmp(ct,'') && strcmp(form,'FAST') % fast ice, which we assume to be 10
    code = 11;
    conc = 10;
    rang = 0;
elseif strcmp(ct,'') && strcmp(form,'STRIPS') % this is for the weird contingency 17Feb14 or thereafter where concentration is '' and form is 'strips'. Not sure what this means
    code = nan;
    conc = nan;
    rang = nan;
elseif strcmp(ct,'')
    code = 0;
    conc = 0;
    rang = 0;
else
    c1 = str2num(ct(1)); % first tenth value
    c2 = str2num(ct(end-1:end)); % second tenth value
    if c2<10 % if second tenth value is not 10
        code = c1*10+c2; % convert it to a number (i.e. 68 for 6-8)
    else
        c2 = 1;
        code = c1*10+c2; % convert it to a number in special case where second value is 10 (i.e. 91 for 9-10)
        c2 = 10; % retain 10 information of second number for concentration average calculation
    end

    conc = mean([c1 c2]); % compute average concentration
    rang = abs(c2-c1)/2; % compute range of concentration

    if isempty(code)
        code = nan;
        conc = nan;
        rang = nan;
    end
end