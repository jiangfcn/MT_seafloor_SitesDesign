function TopoTrackProfile

%pick out the topographical elevation from a DEM grid along a profile
%defined by the coordinates in the InputParameters function.


% OTHER EXTERNAL FUNCTIONS ARE NEEDED:
% -utm2geo, coordinate converting;
% -geo2utm, coordinate converting;
% Documentation for above two functions:
% http://www.uwgb.edu/dutchs/UsefulData/UTMFormulas.HTM
% http://www.linz.govt.nz/geodetic/conversion-coordinates/projection-conversions/transverse-mercator-preliminary-computations/index.aspx

% -twopoints2line, solving a line equation from two known points;
% -point2verticalline, passing a known point to draw a vertical line of a
% known line (ax+by+c=0), and solving this vertical line's equation;
% -lines2point, solving point's coordinate of two intersecting lines;
% -dis2geopoints, calculating distance between two geographical points on
% the earth.
% -writeout, to write a .txt files out


%
% Copyright 2020-2021
%
% Jiang Feng,
% South China Sea Institute of Oceanology, Chinese Academy of Sciences.
% Email: 381908220@qq.com 
%
% License:
%
% It is a free script; you can modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation; either
% version 2 of the License, or any later version at your option. this
% program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. 

% NOTE: GMT has called in this script.It is an open source collection of
% tools and released under the GNU Lesser General Public % License.
% http://gmt.soest.hawaii.edu/projects/gmt

% Please contact the author if you have any problems on the license.


%
% Input and Parameters from FUNCTION: InputParameters.m
%
u = InputParameters; % a input parameters setting by hand.


%-------------------------------------------------------------------------------
% function is below that you do not need any changes!
%-------------------------------------------------------------------------------

%coord,lon,lat
cent_lon = min(u.coord(:,1)) + 0.5*(max(u.coord(:,1))-min(u.coord(:,1)));
cent_lat = min(u.coord(:,2)) + 0.5*(max(u.coord(:,2))-min(u.coord(:,2)));
[coordxy_x,coordxy_y] = geo2utm(u.coord(:,1),u.coord(:,2),cent_lon,cent_lat); %unit: m
coordxy = [coordxy_x,coordxy_y];

nline = size(coordxy,1); % amount of the segments.

%line P
pline = zeros(size(coordxy,1)-1,3);
for ipt = 1:size(coordxy,1)-1
    pline(ipt,:) = twopoints2line(coordxy(ipt,1),coordxy(ipt+1,1),coordxy(ipt,2),coordxy(ipt+1,2));
end

%Line Q,constrain a belt from line P
qline_t1 = pline(:,3)+u.L*sqrt(pline(:,1).^2+pline(:,2).^2);
qline_t2 = pline(:,3)-u.L*sqrt(pline(:,1).^2+pline(:,2).^2);
qline01 = [pline(:,1),pline(:,2),qline_t1];
qline02 = [pline(:,1),pline(:,2),qline_t2];


%calculating a series of sampling points on the line P
% and storing them

dist = 0; %assign 0 for the initial distance；
coordxy_px = []; coordxy_py = [];coordxy_px_c = []; coordxy_py_c = [];
npt = [];%amount of points in each segment.
for ipt = 1:size(coordxy,1)-1
    coordxy_px = [coordxy_px;coordxy(ipt,1)];
    coordxy_py = [coordxy_py;coordxy(ipt,2)];
    coordxy_px_c = [coordxy_px_c;coordxy(ipt,1)]; %中心点也将从线段的端点开始。
    coordxy_py_c = [coordxy_py_c;coordxy(ipt,2)];
    if coordxy(ipt,1)<coordxy(ipt+1,1) %longitude_previous < longitude_after
        if coordxy(ipt,2)<coordxy(ipt+1,2) %the second point fall in the first quadrant of coordinate with the first point as the center
            slope = atan(abs(pline(ipt,1)));
            interv_x = u.interv*cos(slope);
            interv_y = u.interv*sin(slope);
            while coordxy_px(end)+interv_x<coordxy(ipt+1,1)&&coordxy_py(end)+interv_y<coordxy(ipt+1,2)
                coordxy_pxtmp = coordxy_px(end) + interv_x;
                coordxy_pxtmp_c = coordxy_px(end) + interv_x/2;
                coordxy_pytmp = coordxy_py(end) + interv_y;
                coordxy_pytmp_c = coordxy_py(end) + interv_y/2;
                coordxy_px = [coordxy_px;coordxy_pxtmp];
                coordxy_py = [coordxy_py;coordxy_pytmp];
                coordxy_px_c = [coordxy_px_c;coordxy_pxtmp_c];
                coordxy_py_c = [coordxy_py_c;coordxy_pytmp_c];
                dist = [dist;dist(end)+u.interv];
            end
            dist_last = sqrt((coordxy(ipt+1,1)-coordxy_px(end))^2 + (coordxy(ipt+1,2)-coordxy_py(end))^2);
            dist = [dist;dist(end)+dist_last];
        elseif coordxy(ipt,2)>coordxy(ipt+1,2) %the second point fall in the first quadrant of coordinate with the first point as the center
            slope = atan(abs(pline(ipt,1)));
            interv_x = u.interv*cos(slope);
            interv_y = -1*u.interv*sin(slope);
            while coordxy_px(end)+interv_x<coordxy(ipt+1,1)&&coordxy_py(end)+interv_y>coordxy(ipt+1,2)
                coordxy_pxtmp = coordxy_px(end) + interv_x;
                coordxy_pxtmp_c = coordxy_px(end) + interv_x/2;
                coordxy_pytmp = coordxy_py(end) + interv_y;
                coordxy_pytmp_c = coordxy_py(end) + interv_y/2;
                coordxy_px = [coordxy_px;coordxy_pxtmp];
                coordxy_py = [coordxy_py;coordxy_pytmp];
                coordxy_px_c = [coordxy_px_c;coordxy_pxtmp_c];
                coordxy_py_c = [coordxy_py_c;coordxy_pytmp_c];
                dist = [dist;dist(end)+u.interv];
            end
            dist_last = sqrt((coordxy(ipt+1,1) - coordxy_px(end))^2 + (coordxy(ipt+1,2)-coordxy_py(end))^2);
            dist = [dist;dist(end)+dist_last];
        else
            slope = atan(abs(pline(ipt,1)));
            interv_x = u.interv*cos(slope);
            interv_y = -1*u.interv*sin(slope);
            while coordxy_px(end)+interv_x<coordxy(ipt+1,1)
                coordxy_pxtmp = coordxy_px(end) + interv_x;
                coordxy_pxtmp_c = coordxy_px(end) + interv_x/2;
                coordxy_pytmp = coordxy_py(end) + interv_y;
                coordxy_pytmp_c = coordxy_py(end) + interv_y/2;
                coordxy_px = [coordxy_px;coordxy_pxtmp];
                coordxy_py = [coordxy_py;coordxy_pytmp];
                coordxy_px_c = [coordxy_px_c;coordxy_pxtmp_c];
                coordxy_py_c = [coordxy_py_c;coordxy_pytmp_c];
                dist = [dist;dist(end)+u.interv];
            end
            dist_last = abs(coordxy(ipt+1,1) - coordxy_px(end));
            dist = [dist;dist(end)+dist_last];
        end
        
    elseif coordxy(ipt,1)>coordxy(ipt+1,1) %lon_prev > lon_after
        if coordxy(ipt,2)<coordxy(ipt+1,2) %lat_prev < lat_after
            slope = atan(abs(pline(ipt,1)));
            interv_x = -1*u.interv*cos(slope);
            interv_y = u.interv*sin(slope);
            while coordxy_px(end)+interv_x>coordxy(ipt+1,1)&&coordxy_py(end)+interv_y<coordxy(ipt+1,2)
                coordxy_pxtmp = coordxy_px(end) + interv_x;
                coordxy_pxtmp_c = coordxy_px(end) + interv_x/2;
                coordxy_pytmp = coordxy_py(end) + interv_y;
                coordxy_pytmp_c = coordxy_py(end) + interv_y/2;
                coordxy_px = [coordxy_px;coordxy_pxtmp];
                coordxy_py = [coordxy_py;coordxy_pytmp];
                coordxy_px_c = [coordxy_px_c;coordxy_pxtmp_c];
                coordxy_py_c = [coordxy_py_c;coordxy_pytmp_c];
                dist = [dist;dist(end)+u.interv];
            end
            dist_last = sqrt((coordxy(ipt+1,1) - coordxy_px(end))^2 + (coordxy(ipt+1,2)-coordxy_py(end))^2);
            dist = [dist;dist(end)+dist_last];
            
        elseif coordxy(ipt,2)>coordxy(ipt+1,2)
            slope = atan(abs(pline(ipt,1)));
            interv_x = -1*u.interv*cos(slope);
            interv_y = -1*u.interv*sin(slope);
            while coordxy_px(end)+interv_x>coordxy(ipt+1,1)&&coordxy_py(end)+interv_y>coordxy(ipt+1,2)
                coordxy_pxtmp = coordxy_px(end) + interv_x;
                coordxy_pxtmp_c = coordxy_px(end) + interv_x/2;
                coordxy_pytmp = coordxy_py(end) + interv_y;
                coordxy_pytmp_c = coordxy_py(end) + interv_y/2;
                coordxy_px = [coordxy_px;coordxy_pxtmp];
                coordxy_py = [coordxy_py;coordxy_pytmp];
                coordxy_px_c = [coordxy_px_c;coordxy_pxtmp_c];
                coordxy_py_c = [coordxy_py_c;coordxy_pytmp_c];
                dist = [dist;dist(end)+u.interv];
            end
            dist_last = sqrt((coordxy(ipt+1,1) - coordxy_px(end))^2 + (coordxy(ipt+1,2)-coordxy_py(end))^2);
            dist = [dist;dist(end)+dist_last];
        else
            slope = atan(abs(pline(ipt,1)));
            interv_x = -1*u.interv*cos(slope);
            interv_y = u.interv*sin(slope);
            while coordxy_px(end)+interv_x>coordxy(ipt+1,1)
                coordxy_pxtmp = coordxy_px(end) + interv_x;
                coordxy_pxtmp_c = coordxy_px(end) + interv_x/2;
                coordxy_pytmp = coordxy_py(end) + interv_y;
                coordxy_pytmp_c = coordxy_py(end) + interv_y/2;
                coordxy_px = [coordxy_px;coordxy_pxtmp];
                coordxy_py = [coordxy_py;coordxy_pytmp];
                coordxy_px_c = [coordxy_px_c;coordxy_pxtmp_c];
                coordxy_py_c = [coordxy_py_c;coordxy_pytmp_c];
                dist = [dist;dist(end)+u.interv];
            end
            dist_last = abs(coordxy_px(end)-coordxy(ipt+1,1));
            dist = [dist;dist(end)+dist_last];
        end
    else %lon_prev = lon_after
        
        if coordxy(ipt,2)<coordxy(ipt+1,2)
            interv_x = 0;
            interv_y = u.interv;
            while coordxy_py(end)+interv_y<coordxy(ipt+1,2)
                coordxy_pxtmp = coordxy_px(end) + interv_x;
                coordxy_pxtmp_c = coordxy_px(end) + interv_x/2;
                coordxy_pytmp = coordxy_py(end) + interv_y;
                coordxy_pytmp_c = coordxy_py(end) + interv_y/2;
                coordxy_px = [coordxy_px;coordxy_pxtmp];
                coordxy_py = [coordxy_py;coordxy_pytmp];
                coordxy_px_c = [coordxy_px_c;coordxy_pxtmp_c];
                coordxy_py_c = [coordxy_py_c;coordxy_pytmp_c];
                dist = [dist;dist(end)+u.interv];
            end
            dist_last = abs(coordxy_py(end)-coordxy(ipt+1,2));
            dist = [dist;dist(end)+dist_last];
        elseif coordxy(ipt,2)>coordxy(ipt+1,2)
            interv_x = 0;
            interv_y = -1*u.interv;
            while coordxy_py(end)+interv_y>coordxy(ipt+1,2)
                coordxy_pxtmp = coordxy_px(end) + interv_x;
                coordxy_pxtmp_c = coordxy_px(end) + interv_x/2;
                coordxy_pytmp = coordxy_py(end) + interv_y;
                coordxy_pytmp_c = coordxy_py(end) + interv_y/2;
                coordxy_px = [coordxy_px;coordxy_pxtmp];
                coordxy_py = [coordxy_py;coordxy_pytmp];
                coordxy_px_c = [coordxy_px_c;coordxy_pxtmp_c];
                coordxy_py_c = [coordxy_py_c;coordxy_pytmp_c];
                dist = [dist;dist(end)+u.interv];
            end
            dist_last = abs(coordxy_py(end) - coordxy(ipt+1,2));
            dist = [dist;dist(end)+dist_last];
        else
            error('there are two same points in the inputed coord');
        end
    end
    npt = [npt;length(coordxy_px)]; %the amount of sampling points on each line P.
end
coordxy_px = [coordxy_px;coordxy(end,1)]; %add the last point's x
coordxy_py = [coordxy_py;coordxy(end,2)]; %add the last point's y
npt(end) = npt(end)+1; %counting the last point into the end segment's amount of points.

coordxy_px_c = [coordxy_px_c;coordxy(end,1)]; %adding the last point's x as center point's x
coordxy_py_c = [coordxy_py_c;coordxy(end,2)]; %adding the last point's y as center point's y

i=1;
for icross = 1:length(coordxy_px_c)
    if icross <= npt(i)
        Cline(icross,:) = point2verticalline(coordxy_px_c(icross),coordxy_py_c(icross),pline(i,1),pline(i,2));
    else
        i=i+1;
        Cline(icross,:) = point2verticalline(coordxy_px_c(icross),coordxy_py_c(icross),pline(i,1),pline(i,2));
    end
    if i>size(pline)
        error('the maximum of icross is wrong');
    end
end

%calculating intersection points of vertical lines and line Q
j =1;
for inode = 1:length(Cline)
    if inode <= npt(j)
        nodep01(inode,:) = lines2point(qline01(j,1),qline01(j,2),qline01(j,3),Cline(inode,1),Cline(inode,2),Cline(inode,3));
        nodep02(inode,:) = lines2point(qline02(j,1),qline02(j,2),qline02(j,3),Cline(inode,1),Cline(inode,2),Cline(inode,3));
    else
        j=j+1;
        nodep01(inode,:) = lines2point(qline01(j,1),qline01(j,2),qline01(j,3),Cline(inode,1),Cline(inode,2),Cline(inode,3));
        nodep02(inode,:) = lines2point(qline02(j,1),qline02(j,2),qline02(j,3),Cline(inode,1),Cline(inode,2),Cline(inode,3));
    end
    if j > size(pline)
        error('the maximum of inode is wrong');
    end
end


% plot(coordxy_px,coordxy_py,'r.');hold on
% plot(coordxy(:,1),coordxy(:,2),'ks','markersize',20);

figure(1) % showing the sampling points and lines
for j = 1:size(Cline,1)
    xj1 = coordxy_px_c(j)+40000; xj2 = coordxy_px_c(j)-40000;
    yj1 = -1*(Cline(j,1)*xj1+Cline(j,3))/Cline(j,2);
    yj2 = -1*(Cline(j,1)*xj2+Cline(j,3))/Cline(j,2);
    plot([xj1;xj2],[yj1;yj2],'k-');hold on
end
axis equal

% plot(nodep01(:,1),nodep01(:,2),'k.');hold on
% plot(nodep02(:,1),nodep02(:,2),'b.');
% hold off;

%converting the xyz of nodes to geographical coordinates.
[nodeplon01,nodeplat01] = utm2geo(nodep01(:,1),nodep01(:,2),cent_lon,cent_lat);
[nodeplon02,nodeplat02] = utm2geo(nodep02(:,1),nodep02(:,2),cent_lon,cent_lat);
nodeplonlat01 = [nodeplon01,nodeplat01]; nodeplonlat02 = [nodeplon02,nodeplat02];

[coordxy_lon,coordxy_lat] = utm2geo(coordxy_px,coordxy_py,cent_lon,cent_lat);

%estimate an approximated range of DEM we need.
maxdeg = sqrt(2)*u.L/110000;


%% clear memory
clear coordxy_px_c coordxy_py_c coordxy_px coordxy_py ...,
    coordxy_pytmp_c coordxy_pxtmp_c coordxy_pxtmp coordxy_pytmp ...,
    nodeplon01 nodeplon02 nodeplat01 nodeplat02 nodep01 nodep01 ...,
    Cline ...,
    

%% calculating elevation from DEM

DEMv_all = [];
% all the information of the profile, including:
% 'lon, lat, distance, sampling points' elevation, maximum, minimum, and mean elevation
for ipt02 = 1:length(coordxy_lon)-1
    R = [coordxy_lon(ipt02)-maxdeg,coordxy_lon(ipt02)+maxdeg,coordxy_lat(ipt02)-maxdeg,coordxy_lat(ipt02)+maxdeg];  % xmin/xmax/ymin/ymax
    
    %     system(['gmt grdcut ', DEM, ' -Gtmp.grd -R',num2str(R(1)),'/',num2str(R(2)),'/',num2str(R(3)),'/',num2str(R(4))]);
    %     system(['gmt grd2xyz tmp.grd > tmp.xyz']);
    warning('off');
    tmpgrd = gmt(['grdcut  ',u.DEM,' -R',num2str(R(1)),'/',num2str(R(2)),'/',num2str(R(3)),'/',num2str(R(4))]);
    tmpxyz = gmt('grd2xyz  ',tmpgrd); % a structure
    gmt('destroy');
    
    polyg = [nodeplonlat01(ipt02,1),nodeplonlat01(ipt02,2);
        nodeplonlat02(ipt02,1),nodeplonlat02(ipt02,2);
        nodeplonlat02(ipt02+1,1),nodeplonlat02(ipt02+1,2);
        nodeplonlat01(ipt02+1,1),nodeplonlat01(ipt02+1,2);];
    
    tmpxyz = tmpxyz.data;
    [in, on] = inpolygon(tmpxyz(:,1),tmpxyz(:,2),polyg(:,1),polyg(:,2));
    if any(in)
        data_p = [tmpxyz(logical(in+on),1),tmpxyz(logical(in+on),2),tmpxyz(logical(in+on),3)];
    else
        error('Didn''t find any points in a samll topographical belt, please assign a bigger value for u.L or u.interv');
    end
    clear tmpxyz  tmpgrd
    DEMmax = max(data_p(:,3)); DEMmin = min(data_p(:,3)); DEMmean = mean(data_p(:,3));
    
    %interpolation from DEM and calculating elevation.
    [dp,indp] = sort(dis2geopoints(data_p(:,1),data_p(:,2),coordxy_lon(ipt02),coordxy_lat(ipt02)));
    drecip = 1./dp;
    %     dsum = sum(drecip(1:8));
    dsum = sum(drecip(1));
    %DEMp = sum(drecip(1:8).*data_p(indp(1:8),3))/dsum;
    DEMp = sum(drecip(1).*data_p(indp(1),3))/dsum;
    %     lonx = unique(data_p(:,1)); laty = unique(data_p(:,2)); [dmesh_x,dmesh_y] = meshgrid(lonx,laty);
    %     DEMp = interp2(data_p(:,1),data_p(:,2),data_p(:,3),coordxy_lon(ipt),coordxy_lat(ipt));
    DEMv = [DEMp,DEMmax,DEMmin,DEMmean];
    DEMv_all = [DEMv_all;DEMv];
    clear DEMv in on data_p DEMmax DEMmin DEMmean DEMp
    % delete tmp.grd tmp.xyz
end
%% plotting the elevations
figure(2);
h1 = plot(dist(1:end-1)/1000,DEMv_all(:,2)/1000,'r-'); hold on
h2 = plot(dist(1:end-1)/1000,DEMv_all(:,3)/1000,'b-');
h3 = plot(dist(1:end-1)/1000,DEMv_all(:,4)/1000,'k-');

%% adding key points from a text file.
if exist(u.KeyPointsFile,'file')
    if strcmp(u.KeyPointsFile(end-3:end),'.txt') %Key points txt file
        fid = fopen(u.KeyPointsFile);
        linstr = fgetl(fid); i=1;
        while linstr ~=-1
            textlonlat(i,:) = sscanf(linstr,'%f%f');
            strtmp = sscanf(linstr,'%s');
            ptextname{i,:} = strtmp(strfind(strtrim(strtmp),'/')+1:end);
            i = i+1;
            linstr = fgetl(fid);
        end
    elseif strcmp(u.KeyPointsFile(end-3:end),'.kml') %Key points kml file
        textstruct = gmt(['kml2gmt ',u.KeyPointsFile,' -Fs -Vq']);
        textlonlat = [];
        fid_0 = fopen([u.KeyPointsFile(1:end-4),'FromKML.txt'],'wt');
        header_0 = '#ID  #lon  #lat';
        fprintf(fid_0,'%s\n',header_0);
        for ipt = 1:length(textstruct)
            header = textstruct(ipt).header;
            ind = strfind(header,'"');
            ptextname{ipt,:} = header(ind(1)+1:ind(2)-1);
            tmp = textstruct(ipt).data;
            textlonlat = [textlonlat;tmp];
            fprintf(fid_0,'%s\t',ptextname{ipt,:});
            fprintf(fid_0,'%16.8f  %16.8f\n',textlonlat(ipt,:));            
        end
        fclose(fid_0);

    else
        error('!!!the Key Points File should be .txt or .kml!!!');
    end
    
    %Calculating points' distance and elevation after project on the line P!
    % so, calculation below is not the real elevation for the key points in the text file, it is
    % a elevation of point on the line P which is the closest one to the key point!
    [textx,texty] = geo2utm(textlonlat(:,1),textlonlat(:,2),cent_lon,cent_lat);
    textxy = [textx,texty];
    for icontrl = 1:size(textlonlat,1)
        textdist = dis2geopoints(textlonlat(icontrl,1),textlonlat(icontrl,2),coordxy_lon,coordxy_lat);
        [vt,indt] = min(textdist);
        [sorttmp,indsort] = sort(unique([npt;indt]));
        indt_loc = find(indsort == length(npt)+1); %Choosing a line for the point.
        indt_loc = 1;
        textline = point2verticalline(textxy(icontrl,1),textxy(icontrl,2),pline(indt_loc,1),pline(indt_loc,2));
        ptextxy = lines2point(textline(1),textline(2),textline(3),pline(indt_loc,1),pline(indt_loc,2),pline(indt_loc,3));
        [ptextlon,ptextlat] = utm2geo(ptextxy(1),ptextxy(2),cent_lon,cent_lat);% lon, lat
        
        ptextd = dis2geopoints(ptextlon,ptextlat,coordxy_lon,coordxy_lat);
        [ptextds,ptextd_ind] = sort(ptextd);
        dist_v = min(dist(ptextd_ind(1)),dist(ptextd_ind(2))); %choosing the nearest two points.
        if dist_v == dist(ptextd_ind(1));  ind_v=1;       else  ind_v=2;      end
        ptext_xy(icontrl,1) = dist_v + dis2geopoints(ptextlon,ptextlat,coordxy_lon(ptextd_ind(ind_v)),coordxy_lat(ptextd_ind(ind_v)) );
        % profile distance = left hand points' distance + the length between
        % them, this profile distance is good calculation for sites along a
        % profile!!!!
        
        ptext_xy(icontrl,2) = min([DEMv_all(ptextd_ind(1),2), DEMv_all(ptextd_ind(2),2)]); % choose mean for water depth!!!
        % calculating a mean value of the two mean elevation (4) among
        % topography belt and assigned as the key point's elevation.
    end
    
    fid_1 = fopen([u.KeyPointsFile(1:end-4),'_KeyPointInfo.txt'],'wt');
    header_1 = '#ID  #lon  #lat  #dist  #elev';
    fprintf(fid_1,'%s\n',header_1);
    for iplot = 1:length(ptext_xy)
        plot(ptext_xy(iplot,1)/1000,ptext_xy(iplot,2)/1000,' vk','markeredgecolor','k','markerfacecolor','k','markersize',u.markersize);hold on
        if ~isempty(ptextname)
            text(ptext_xy(iplot,1)/1000 + u.shift_right,ptext_xy(iplot,2)/1000+u.shift_upward,ptextname{iplot}(u.textind:end),...
                'verticalAlignment','top','horizontalAlignment','left','Rotation',90,'FontWeight','bold','Fontsize',u.fontsize_text);
        end
        fprintf(fid_1,'%s\t',ptextname{iplot});
        fprintf(fid_1,'%16.8f  %16.8f  %10.2f  %10.0f\n',[textlonlat(iplot,:),ptext_xy(iplot,:)]);
    end
    ylabel('Elevation(km)');xlabel('Distance(km)');
    set(gca,'fontsize',u.fontsize_gca);
    grid on;
    legend([h1,h2,h3],'max','min','ave'); hold off
    fclose(fid_1);
else
    disp(['Didnot find any KeyPoints in this run!']);
end % adding Key Points File.
hold off; % for figure 2

header_2 = '#lon lat dist p_elev max_elev min_elev mean_elev';
fid_2 = fopen('DEMtrack_profile_info.txt','wt'); writeout(fid_2,[coordxy_lon(1:end-1),coordxy_lat(1:end-1),dist(1:end-1),DEMv_all],header_2);


delete gmt.history

end
