close all; clear all; clc;
% 
% generating equal spacing points on a geospheric line. 
%
n = 20; % number of points you want.
prefix = 'LYT';% provide sites' name
snumb = 2300 + [1:n];
P1 = [15.230000, 107.104293]; % the two points on the endd of a line.
P2 = [17.240000, 107.748574];


%% NOT CAHNGE ANYTHING BELOW=======================================
filename = [prefix,num2str(n),'.txt']; % files will be output
if exist(filename,'file')
    error('The File Already Exist in This Folder');
end
for iid = 1:length(snumb)
    ID(iid) = string([prefix,num2str(snumb(iid))]); 
end

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


