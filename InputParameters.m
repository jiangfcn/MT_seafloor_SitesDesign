function [u] = InputParameters
% Parameters input to TopoTrackProfile.m.

% DEM file and its absolute path, the DEM file should be read by the GMT
% functions.
u.DEM = 'A:\InterpretData\DEM\ETOPO1\ETOPO1_Bed.grd';

% A file contains key points which are needed to be plotted on the DEM
% track profile. Comment the u.KeyPointsFile line if there is no this
% file.
% This file can be a .txt or a .kml from the Google Earth.
% for the .txt: first and second column - longitude and latitude,following
% optional can have the points' ID. It is looks like:
% longitude01  latitude01  / ID01
% longitude02  latitude02  / ID02
% longitude03  latitude03  / ID03
%            .
%            .
%            .
u.KeyPointsFile = 'OBS2016-2POBSPlan_2023.kml';

% Geographical coordinates in two columns to define the profiles. These profiles can form polyline by connecting end-to-end, 
% 1-longitude and 2-latitude (degree!),
u.coord = [140.9445303218392,12.72727761332462;...
           143.003458874507,10.37698870474865]; 

% parallel distance between the center line P and the offset line Q, the topographical belt width will be 2*u.L.
% these two values should be bigger than the resolution of the DEM data.
u.L = 2000; % unit:m
u.interv = 2000; % unit: m

%text the Key points's name on each symbol
u.markersize = 10;
u.shift_right = -5; % unit:km
u.shift_upward = 0.1; % unit:km
u.fontsize_text = 16; % fontsize for the text label of key point.
u.fontsize_gca = 16; % fontsize in the figure of topographical profile
u.textind = 1; % a number from which the key points' name will be shown.
end
