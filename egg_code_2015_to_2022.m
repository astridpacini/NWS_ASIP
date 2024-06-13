function [code,conc,rang] = egg_code_2015_to_2022(ct,form,eggcol)
%PARSE egg code here. This uses the onfomration in NSIDC usergguide to
%assign concentration intervals to an egg code value, and then convert that
%to a sea ice concentration (in tenths). 
%Note that this function handles the post Oct 1 2015 NWS ASIP datafiles, 
%which tabulate concentration information with a character string without a
%dash (e.g. '68', '91'). 
%Written by A. Pacini, apacini@uw.edu, May 2024

if strcmp(ct,'') 
    code = 0;
    conc = 0;
    rang = 0;
elseif strcmp(form,'08') || strcmp(eggcol,'93')% note this takes presedence over ice 
    % concentration information (as is the case in NWS figures)
    code = 11;
    c1 = str2num(ct(1));
    c2 = str2num(ct(2));

    if c1==9 && c2==2 % cases where code is 92, 10/10 ice no range
        conc = 10;
        rang = 0;
    elseif c2 == 1 && c1>0 % cases where code has 10 in second value (e.g. 91)
        rang = abs(10-c1)/2;
        conc = mean([c1 10]);
    else % cases that do not have maximum ice (i.e. 1 or 2 in second digit) but are still fast ice (as denoted by form or eggcol above)
        rang = abs(c2-c1)/2;
        conc = mean([c1 c2]);
    end
else
    code = str2num(ct);
    c1 = str2num(ct(1));
    c2 = str2num(ct(2));
  
    if  code == 92 % cases with heavy ice, not labeld fast ice
        conc = 10;
        rang = 0;
    elseif c2 == 1 && c1>0 % in scenario where second digit indicates 10 (e.g. 91)
        rang = abs(10-c1)/2; % compute range of concentration
        conc = mean([c1 10]); % compute average concentration
    else % in scenario where second digit does not indicate 10
        rang = abs(c2 - c1)/2; % compute range of concentration
        conc = mean([c1 c2]); % compute average concentration
    end
    
end