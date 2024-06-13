%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes in raw National Weather Service Alaska Sea Ice Program
% data and parses, projects, and grids it. The details are provided below,
% but shapefiles are injested, and .mat and .nc files are output.
%
% This script does the following:
% 1. define a gridding domain and resolution.
% 2. read and project the daily shapefiles onto lat/lon domain
% 3. iterate through each polygon in each day and find the grid cells that
% this polygon covers
% 4. fill out the grid information with the polygon data (SIC, area of
% polygon).
% 5. Deal with a few scenarios:
% 6. output gridded fiels containing the following variables:
%       a. LAT
%       b. LON
%       c. date (matlab datetime format)
%       d. iceconc (tenths): sea ice concentration
%       e. rangeconc (tenths): range of sea ice concentration, as defined
%       by polygon range (standard egg code)
%       f. eggcode (#): origical egg code number defining SIC in polygon,
%       see WMO literature (e.g. WMO, 1970)
%       g. flag (#): this is an accounting variable that identifies
%       instances where one polygon was embedded inside another polygon,
%       smaller polygon is always trusted.
%       h. amt_poly (#): the number of polygons in a given day
%
% A few important notes about conventions:
% 1. In 2015, NWS ASIP changed its file convention. Starting Oct. 1, 2015,
% ASIP moved to daily subdirectories, instead of monthly subdirectories.
% 2. ASIP provided 3x per wee maps from March 2007-July 2014 (M/W/Fri), and
% daily maps after that.
% 3. CHECK ON EGG CODE CONVENTION CHANGE--2014 or 2015?
%
% Dependencies: convert the combination of egg code CT (SIC) value and
% form information to a single value for sea ice concentraiton
% egg_code_2007_to_2015: conversion for data prior to Oct. 1, 2015
% egg_code_2015_to_2022: conversion for data post Oct. 1, 2015
%
% Written by A. Pacini, apacini@uw.edu
% Date Created: 4/26/2023
% Date Modified: 5/24/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;


% define grid (this case, 0.05 degree spacing in lat and lon)
lats = [55.125:0.05:80.1250]; % latitude grid
lons1 = [174.75:0.05:179.96]; % longitudes to the west of dateline
lons2 = [-180:0.05:-135.125]; % longitudes to the east of dateline
lons = [lons1,lons2]; % longitude grid (to merge + values west of dateline, - values east of dateline)

[LON,LAT]=meshgrid(lons,lats); % create 2D gridded field
LON = wrapTo360(LON); % use 360 degree convention because dealing with dateline

% year(s) that script should grid
yrs = [2015];

for i = 1:length(yrs) % for each year
    cd(['//Users/astridpacini/Documents/UW_postdoc/research/data/sic/nws/' num2str(yrs(i))]);

    % Get a list of all files and folders in this folder. This is to
    % identify available months
    files = dir;

    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];

    % Extract only those that are directories.
    subFolders = files(dirFlags); % structure with extra info.

    % Get only the folder names into a cell array.
    subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..

    % daily vector for year
    date = datenum([datetime(yrs(i),01,01,0,30,0):caldays(1):datetime(yrs(i),12,31,0,30,0)]); % create daily datetime vector, valid at 00:30

    % pre-populate NaN matrices for each variable for the year
    iceconc = nan(size(LAT,1),size(LAT,2),length(date)); % 3D (lats = rows, lons = columns, time = 3rd dimension)
    rangeconc = iceconc; % 3D
    eggcode = iceconc; % 3D
    flag = iceconc; % 3D
    area = iceconc; % 3D
    amt_poly = nan(1,length(date)); % 2D, one value per day

    % initiate dateT
    dateT = datenum(yrs(i),1,01);

    for j = 1:length(subFolderNames)
        cd(['//Users/astridpacini/Documents/UW_postdoc/research/data/sic/nws/' num2str(yrs(i)),'/' subFolderNames{j}]);

        % now it makes a difference which year we're in. For 2007-2015, we have monthly subdirectories, for 2016-2022, we have daily subdirectories

        % monthly subdirectories until 10/1/2015
        if dateT < datenum(2015,09,30) % last date before the subdirectory switch (this corresponds to the last loop's date)

            % find the shapefiles
            dys = dir('*.shp');

            % now extract information on which days have ice data
            for k = 1:length(dys)
                day = dys(k).name % let this print out, to know at what point the code is

                dateT = datenum(yrs(i),str2num(subFolderNames{j}),str2num(day(end-10:end-9))); % extract date information from filename

                % match to the date of original matrix, thus leaving others nan-ed out
                [~,id] = min(abs(dateT-date));

                % now let's read in data and grid it for this date
                file = day(1:end-4);
                info = shapeinfo(file);
                S = shaperead(file);
                proj = info.CoordinateReferenceSystem; % load projection information from file

                % initialize concentration and flag map
                icec = nan(size(LAT,1),size(LAT,2));
                code = icec;
                rang = icec;
                areaT = icec;
                flagT = icec;

                % now iterate through each polygon
                for ii = 1:length(S)

                    ct_egg = (S(ii).CONCENTRAT); % ice concentration
                    form = S(ii).FORM; % primary form of ice
                    if dateT < datenum(2011,10,14) % capitalization changed convention after 2011
                        area_S = S(ii).SHAPE_Area; % area of polygon
                    else
                        area_S = S(ii).Shape_Area; % area of polygon
                    end

                    x = S(ii).X;
                    y = S(ii).Y;

                    [lat,lon] = projinv(proj,x,y); % this function converts
                    % x and y coordiantes to lat and lon values specified
                    % by the projection information, which is a projcrs
                    % object identified above

                    lon = wrapTo360(lon); % wrap longitude to 360 deg to match grid

                    % rules for egg code convention change:
                    % On Oct 1, 2015 we move from interval into code form,
                    % i.e. the shapefile simply contains egg code assignment
                    % (e.g. 71,78, etc.) instead of interval
                    % ('7-10, 7-8', etc.) (thus switch to 2015-2022
                    % function)

                    if dateT>=datenum(2015,10,01)
                        eggcol = S(ii).COLOR;
                        [ce,ct,rg] = egg_code_2015_to_2022(ct_egg,form,eggcol);
                    else
                        [ce,ct,rg] = egg_code_2007_to_2015(ct_egg,form);
                    end

                    % now for the parsing/separating of polygons
                    idx = find(isnan(lat)) ;     % find positions of NaNs
                    idx = [1 idx length(lat)] ;  % append first and last position, IS THIS RIGHT?

                    % go through each polygon
                    for jj = 1:length(idx)-1
                        pos = idx(jj):idx(jj+1) ;  % get the required position
                        xi = lon(pos);
                        yi = lat(pos);  % get the corodinates

                        % Remove NaN's
                        xi(isnan(xi)) = [] ;
                        yi(isnan(yi)) = [] ;

                        %xi = wrapTo360(xi); % wrap longitude to 360 deg to match grid

                        % if we have a polygon with ice information
                        if ~isempty(ct)
                            in = inpolygon(LAT,LON,yi,xi); % create a logical matrix that tells us where
                            % polygon falls on the grid (1 = gridcells within polygon, 0 = gridcells outside of polygon)
                            inct = in.*ct; % assign each gridcell point in the polygon with SIC value
                            ince = in.*ce; % assign each gridcell point in the polygon with eggcode value
                            inrg = in.*rg; % assign each gridcell point in the polygon with range value
                            inar = in.*area_S; % assign each gridcell point in the polygon with the area of that polygon

                            % find everything outside polygon
                            ou = in == 0;

                            % and make it nan (instead of 0)
                            inct(ou) = nan;
                            ince(ou) = nan;
                            inrg(ou) = nan;
                            inar(ou) = nan;

                            % identify gridcells with new information
                            indnew = inct>0;

                            % identify original value of those gridcells
                            indold = icec>0 & in;

                            % do we have overlapping gridcells?
                            if sum(sum(indnew + indold == 2))>0

                                indtot = indnew+indold;
                                index_overlap = indtot == 2; % the overlapping gridcells

                                % identify gridcells corresponding to old
                                % polygon
                                oldvalue = icec(in);
                                oldarea = areaT(in);

                                % now check to make sure we don't have multiple
                                % values (i.e. multiple polygon overlaps)
                                [unq,IA,~] = unique(oldarea(~isnan(oldarea)));
                                oldvalue = oldvalue(~isnan(oldarea));
                                unqv = oldvalue(IA); % SIC in places with unique area values within polygon

                                % ok, now we need to decide which polygon
                                % is bigger, the old one or the new one. We
                                % will always trust the smaller polygon
                                % (which is generally the polygon loaded in
                                % later). How does this play out?

                                % for every gridcell of overlap, ask what
                                % is larger, the old or the new polygon?

                                % if new polygon is larger than old
                                % polygon, make new poygon have nan and
                                % keep old value

                                % if new polygon is smaller than old
                                % polygon, make old polygon have nan and
                                % apply new value

                                num_overlap = sum(sum(index_overlap)); % the amount of overlapping gridcells
                                index_overlap_matrix = find(index_overlap); % the actual coordinates in the matrix of overlap

                                for ijk = 1:num_overlap
                                    oldarea = areaT(index_overlap_matrix(ijk));
                                    newarea = inar(index_overlap_matrix(ijk));
                                    if newarea<oldarea % take the new info (make old info nan)
                                        icec(index_overlap_matrix(ijk)) = nan;
                                        code(index_overlap_matrix(ijk)) = nan;
                                        rang(index_overlap_matrix(ijk)) = nan;
                                        areaT(index_overlap_matrix(ijk)) = nan;
                                    elseif newarea>oldarea % take the old info (make new info nan)
                                        inct(index_overlap_matrix(ijk)) = nan;
                                        ince(index_overlap_matrix(ijk)) = nan;
                                        inrg(index_overlap_matrix(ijk)) = nan;
                                        inar(index_overlap_matrix(ijk)) = nan;
                                    elseif newarea==oldarea % take the new info (make old info nan)
                                        icec(index_overlap_matrix(ijk)) = nan;
                                        code(index_overlap_matrix(ijk)) = nan;
                                        rang(index_overlap_matrix(ijk)) = nan;
                                        areaT(index_overlap_matrix(ijk)) = nan;
                                    end
                                    flagT(index_overlap_matrix(ijk)) = sum(flagT(index_overlap_matrix(ijk)),1,'omitnan'); % keep tabs on where we saw overlap
                                end
                            end

                            % now merge the polygon into the broader
                            % gridded fields

                            % make a 3D vector, with dimensions being LAT,
                            % LON, old grid for that day, updated grid for
                            % that day, with only that polygon's worth of
                            % info

                            % SIC
                            c = cat(3,icec,inct);

                            nan_check1c = isnan(icec); % find old nans
                            nan_check2c = isnan(inct); % find new nans
                            nan_check_c = (nan_check1c+nan_check2c)==2; % find where both are nan

                            icec = sum(c,3,'omitnan'); % turn 3D to 2D matrix
                            icec(nan_check_c) = nan; % re-insert nans where we have no information

                            % egg code value
                            d = cat(3,code,ince);

                            nan_check1d = isnan(code); % find old nans
                            nan_check2d = isnan(ince); % find new nans
                            nan_check_d = (nan_check1d+nan_check2d)==2; % find where both are nan

                            code = sum(d,3,'omitnan'); % turn 3D to 2D matrix
                            code(nan_check_d) = nan; % re-insert nans where we have no information

                            % SIC range
                            r = cat(3,rang,inrg);

                            nan_check1r = isnan(rang); % find old nans
                            nan_check2r = isnan(inrg); % find new nans
                            nan_check_r = (nan_check1r+nan_check2r)==2; % find where both are nan

                            rang = sum(r,3,'omitnan'); % turn 3D to 2D matrix
                            rang(nan_check_r) = nan; % re-insert nans where we have no information

                            % area
                            a = cat(3,areaT,inar);

                            nan_check1a = isnan(areaT); % find old nans
                            nan_check2a = isnan(inar); % find new nans
                            nan_check_a = (nan_check1a+nan_check2a)==2; % find where both are nan

                            areaT = sum(a,3,'omitnan'); % turn 3D to 2D matrix
                            areaT(nan_check_a) = nan; % re-insert nans where we have no information
                        end

                    end
                end

                      % now save the information for that day
                iceconc(:,:,id) = icec;
                rangeconc(:,:,id) = rang;
                eggcode(:,:,id) = code;
                flag(:,:,id) = flagT;
                area(:,:,id) = areaT;
                amt_poly(id) = length(S);
                clear icec rang code flagT areaT S

            end
        else % daily subdirectories

            % identify available days
            files = dir;

            % Get a logical vector that tells which is a directory.
            dirFlags = [files.isdir];

            % Extract only those that are directories.
            subsubFolders = files(dirFlags); % structure with extra info.

            % Get only the folder names into a cell array.
            subsubFolderNames = {subsubFolders(3:end).name}; % Start at 3 to skip . and ..

            % now cycle
            for k = 1:length(subsubFolderNames)
                cd(['//Users/astridpacini/Documents/UW_postdoc/research/data/sic/nws/' num2str(yrs(i)),'/' subFolderNames{j},'/' subsubFolderNames{k}]);
                dayinfo = subsubFolderNames{k} % let this print out, to know at what point the code is

                dateT = datenum(yrs(i),str2num(subFolderNames{j}),str2num(dayinfo(end-1:end))); % extract date information from filename

                % match to the date of original matrix, thus leaving others
                % nan-ed out
                [~,id] = min(abs(dateT-date));

                % now let's read in data and grid it for this date
                file = dayinfo;
                info = shapeinfo(file);
                S = shaperead(file);
                proj = info.CoordinateReferenceSystem; % load projection information from file

                % initialize concentration and flag map
                icec = nan(size(LAT,1),size(LAT,2));
                code = icec;
                rang = icec;
                areaT = icec;
                flagT = nan(size(LAT,1),size(LAT,2));

                % now iterate through each polygon
                for ii = 1:length(S)
                    ct_egg = (S(ii).CT); % ice concentration
                    form = S(ii).FP; % primary form of ice
                    eggcol = S(ii).COLOR;
                    area_S = S(ii).Shape_Area;

                    x = S(ii).X;
                    y = S(ii).Y;

                    [lat,lon] = projinv(proj,x,y); % this function converts
                    % x and y coordiantes to lat and lon values specified
                    % by the projection information, which is a projcrs
                    % object identified above

                    lon = wrapTo360(lon); % wrap longitude to 360 deg to match grid

                    % rules for egg code
                    [ce,ct,rg] = egg_code_2015_to_2022(ct_egg,form,eggcol); % this outputs code

                    % now for the parsing/separating of polygons
                    idx = find(isnan(lat)) ;     % find positions of NaNs
                    idx = [1 idx length(lat)] ;  % append first and last position

                    % go through each polygon
                    for jj = 1:length(idx)-1
                        pos = idx(jj):idx(jj+1) ;  % get the required position
                        xi = lon(pos);
                        yi = lat(pos);  % get the corodinates

                        % Remove NaN's
                        xi(isnan(xi)) = [] ;
                        yi(isnan(yi)) = [] ;

                        % if we have a polygon with ice information
                        if ~isempty(ct)
                            in = inpolygon(LAT,LON,yi,xi); % create a logical matrix that tells us where
                            % polygon falls on the grid (1 = gridcells within polygon, 0 = gridcells outside of polygon)

                            inct = in.*ct; % assign each gridcell point in the polygon with SIC value
                            ince = in.*ce; % assign each gridcell point in the polygon with egg code value
                            inrg = in.*rg; % assign each gridcell point in the polygon with range value
                            inar = in.*area_S; % assign each gridcell point in the polygon with the area of that polygon

                            % find everything outside polygon
                            ou = in == 0;

                            % and make it nan (instead of 0)
                            inct(ou) = nan;
                            ince(ou) = nan;
                            inrg(ou) = nan;
                            inar(ou) = nan;

                            % identify gridcells with new information
                            indnew = inct>0;

                            % identify original value of those gridcells
                            indold = icec>0 & in;

                            % do we have overlapping gridcells?
                            if sum(sum(indnew + indold == 2))>0

                                indtot = indnew+indold;
                                index_overlap = indtot == 2; % the overlapping gridcells

                                % identify gridcells corresponding to old
                                % polygon
                                oldvalue = icec(in);
                                oldarea = areaT(in);

                                % now check to make sure we don't have multiple
                                % values (i.e. multiple polygon overlaps)
                                [unq,IA,~] = unique(oldarea(~isnan(oldarea)));
                                oldvalue = oldvalue(~isnan(oldarea));
                                unqv = oldvalue(IA); % SIC in places with unique area values within polygon

                                % ok, now we need to decide which polygon
                                % is bigger, the old one or the new one. We
                                % will always trust the smaller polygon
                                % (which is generally the polygon loaded in
                                % later). How does this play out?

                                % for every gridcell of overlap, ask what
                                % is larger, the old or the new polygon?

                                % if new polygon is larger than old
                                % polygon, make new poygon have nan and
                                % keep old value

                                % if new polygon is smaller than old
                                % polygon, make old polygon have nan and
                                % apply new value

                                num_overlap = sum(sum(index_overlap)); % the amount of overlapping gridcells
                                index_overlap_matrix = find(index_overlap); % the actual coordinates in the matrix of overlap

                                for ijk = 1:num_overlap
                                    oldarea = areaT(index_overlap_matrix(ijk));
                                    newarea = inar(index_overlap_matrix(ijk));
                                    if newarea<oldarea % take the new info (make old info nan)
                                        icec(index_overlap_matrix(ijk)) = nan;
                                        code(index_overlap_matrix(ijk)) = nan;
                                        rang(index_overlap_matrix(ijk)) = nan;
                                        areaT(index_overlap_matrix(ijk)) = nan;
                                    elseif newarea>oldarea % take the old info (make new info nan)
                                        inct(index_overlap_matrix(ijk)) = nan;
                                        ince(index_overlap_matrix(ijk)) = nan;
                                        inrg(index_overlap_matrix(ijk)) = nan;
                                        inar(index_overlap_matrix(ijk)) = nan;
                                    elseif newarea==oldarea % take the new info (make old info nan)
                                        icec(index_overlap_matrix(ijk)) = nan;
                                        code(index_overlap_matrix(ijk)) = nan;
                                        rang(index_overlap_matrix(ijk)) = nan;
                                        areaT(index_overlap_matrix(ijk)) = nan;
                                    end
                                    flagT(index_overlap_matrix(ijk)) = sum(flagT(index_overlap_matrix(ijk)),1,'omitnan'); % keep tabs on where we see overlap
                                end
                            end

                            % now merge the polygon into the broader
                            % gridded fields

                            % make a 3D vector, with dimensions being LAT,
                            % LON, old grid for that day, updated grid for
                            % that day, with only that polygon's worth of
                            % info

                            % SIC
                            c = cat(3,icec,inct);

                            nan_check1c = isnan(icec); % find old nans
                            nan_check2c = isnan(inct); % find new nans
                            nan_check_c = (nan_check1c+nan_check2c)==2; % find where both are nan

                            icec = sum(c,3,'omitnan'); % turn 3D to 2D matrix
                            icec(nan_check_c) = nan; % re-insert nans where we have no information

                            % egg code value
                            d = cat(3,code,ince);

                            nan_check1d = isnan(code); % find old nans
                            nan_check2d = isnan(ince); % find new nans
                            nan_check_d = (nan_check1d+nan_check2d)==2; % find where both are nan

                            code = sum(d,3,'omitnan'); % turn 3D to 2D matrix
                            code(nan_check_d) = nan; % re-insert nans where we have no information

                            % SIC range
                            r = cat(3,rang,inrg);

                            nan_check1r = isnan(rang); % find old nans
                            nan_check2r = isnan(inrg); % find new nans
                            nan_check_r = (nan_check1r+nan_check2r)==2; % find where both are nan

                            rang = sum(r,3,'omitnan'); % turn 3D to 2D matrix
                            rang(nan_check_r) = nan; % re-insert nans where we have no information

                            % area
                            a = cat(3,areaT,inar);

                            nan_check1a = isnan(areaT); % find old nans
                            nan_check2a = isnan(inar); % find new nans
                            nan_check_a = (nan_check1a+nan_check2a)==2; % find where both are nan

                            areaT = sum(a,3,'omitnan'); % turn 3D to 2D matrix
                            areaT(nan_check_a) = nan; % re-insert nans where we have no information
                        end
                    end
                end

                % now save the information for that day
                iceconc(:,:,id) = icec;
                rangeconc(:,:,id) = rang;
                eggcode(:,:,id) = code;
                flag(:,:,id) = flagT;
                area(:,:,id) = areaT;
                amt_poly(id) = length(S);
                clear icec rang code flagT areaT S
            end
        end
    end

    % save the data to a mat file
    cd('//users/astridpacini/Documents/UW_postdoc/research/data/sic/nws/gridded_product/mat_files/')
    save(['nws_' num2str(yrs(i))],'LAT','LON','date','iceconc','rangeconc','eggcode','flag','area','amt_poly')

    % save the data to a NetCDF file
    cd('//users/astridpacini/Documents/UW_postdoc/research/data/sic/nws/gridded_product/netCDF_files/')
    [X,Y,Z] = size(iceconc); % define the dimensions

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'LAT','Dimensions',{'x',X,'y',Y},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'LAT',LAT);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'LAT', 'long_name','degrees north');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'LAT', 'units','degrees_north');

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'LON','Dimensions',{'x',X,'y',Y},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'LON',LON);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'LON', 'long_name','degrees');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'LON', 'units','degrees');

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'date','Dimensions',{'z',Z},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'date',date);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'date', 'long_name','time_value');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'date', 'units','matlab_datetime');

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'iceconc','Dimensions',{'x',X,'y',Y,'z',Z},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'iceconc',date);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'iceconc', 'long_name','ice concentration');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'iceconc', 'units','tenths');

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'rangeconc','Dimensions',{'x',X,'y',Y,'z',Z},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'rangeconc',date);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'rangeconc', 'long_name','ice concentration range');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'rangeconc', 'units','tenths');

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'eggcode','Dimensions',{'x',X,'y',Y,'z',Z},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'eggcode',date);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'eggcode', 'long_name','original_egg_value');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'eggcode', 'units','egg_code');

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'flag','Dimensions',{'x',X,'y',Y,'z',Z},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'flag',date);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'flag', 'long_name','overlapping_polygons');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'flag', 'units','flag_for_overlapping_polygons');

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'area','Dimensions',{'x',X,'y',Y,'z',Z},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'area',date);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'area', 'long_name','source_polygon_area');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'area', 'units','?');

    nccreate(['nws_' num2str(yrs(i)) '.nc'],'amt_poly','Dimensions',{'z',Z},'FillValue','disable');
    ncwrite(['nws_' num2str(yrs(i)) '.nc'],'amt_poly',amt_poly);
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'amt_poly', 'long_name','number_of_polygons');
    ncwriteatt(['nws_' num2str(yrs(i)) '.nc'], 'amt_poly', 'units','number');

end