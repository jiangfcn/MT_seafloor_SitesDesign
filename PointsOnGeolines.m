close all; clear all; clc;
% generating points on a geospheric line and making the points distribute with equal spacing. 
n = 5; % number of points you want.
filename = 'A:\Project\Research\MarianaSubz\Doc\MTsites\MR2023locs_cp.txt'; % files will be output
% provide sites' name
% ID = ["MR2301","MR2302","MR2303","MR2304","MR2305",...
%         "MR2306","MR2307","MR2308","MR2309","MR2310",...
%         "MR2311","MR2312","MR2313","MR2314","MR2315"];
%     ,...
ID = ["MR2316","MR2317","MR2318","MR2319","MR2320"]; % sites's name
% P1 = [1.726530424612083,148.2032171070533];
% P2 = [1.92349706825732,150.4759697772407];
% P1 = [12.669000, 140.996111];
% P2 = [11.643547, 141.898852];
P1 = [11.183519, 142.302020]; % the two points on the endd of a line.
P2 = [10.594834, 142.814680];

%% NOT CAHNGE ANYTHING BELOW=======================================
geocent = [min(P1(1),P2(1)) + abs(P1(1)-P2(1)), min(P1(2),P2(2)) + abs(P1(2)-P2(2))];
long = [P1(2);P2(2)];
lat = [P1(1);P2(1)];
cent_long = geocent(2);
cent_lat = geocent(1);
[Px,Py] = geo2utm(long,lat,cent_long,cent_lat);
linepar = twopoints2line(Px(1),Px(2),Py(1),Py(2));
x_row = linspace(min(Px(1),Px(2)),max(Px(1),Px(2)),n);
y_row = linepar(1).*x_row + linepar(3);
[new_long,new_lat] = utm2geo(x_row,y_row,cent_long,cent_lat);

fid = fopen(filename,'wt');
header = '#long   lat   ID';
data = [new_long',new_lat'];
writeout(fid,data,header,ID);


