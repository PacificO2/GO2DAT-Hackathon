% This tool allows the user to input two sets of coordinates with regions
% of interest to contrast against coastal data

% It uses the PacificO2 database, compiled with PacificO2_main.m

% History:
% Created the 3.May.2022 by the PacificO2 team for the GO2Dat Hackaton
% 2022.

%clear all; close all; clc

%% === User defined ===

% Example for Stratus mooring vs Tongoy mooring
%Stratus = [-20,-85]; % Lat, lon of open ocean box
%Tongoy = [-30.275, -71.562]; 

b1_lower_lon = -86;
b1_higher_lon = -84;
b1_lower_lat = -21;
b1_higher_lat = -19;

b2_lower_lon = -72.562;
b2_higher_lon = -70.562;
b2_lower_lat = -31.275;
b2_higher_lat = -29.275;

%% === Load the data ===

d = load('/Users/ana/Data/hackathon_go2dat_2022/Data/Processed/PacificO2_wod_mirai_coastal.mat');

lon = d.d.lon; 
lat = d.d.lat; 

%% === Select the data within the two boxes ===

b1_idx = find(lat >= b1_lower_lat & lat <= b1_higher_lat & ...
    lon >= b1_lower_lon & lon <= b1_higher_lon);

b2_idx = find(lat >= b2_lower_lat & lat <= b2_higher_lat & ...
    lon >= b2_lower_lon & lon <= b2_higher_lon);

%% === Make a figure with the boxes ===

m_proj('mercator','longitudes',[-90 -70],'latitudes',[-35 -15]);
m_gshhs_h('save','gumby');
    
f1 = figure;
m_grid('tickdir','out')
hold on
h = m_scatter(lon,lat,'k.'); % Plot all the PacificO2 data
hold on
h = m_scatter(lon(b1_idx),lat(b1_idx),'bo'); % Plot data from box1
hold on
h = m_scatter(lon(b2_idx),lat(b2_idx),'ro'); % Plot the Coastal moorings
hold on
m_usercoast('gumby','patch',[0.8 0.8 0.8],'edgecolor','none');


% === Save the figure ===
figdir = '/Users/ana/Data/hackathon_go2dat_2022/figures/';
if ~exist(figdir, 'dir')
    mkdir(figdir)
end

filename1 = ['selected_boxes'];
% for eps
set(gcf,'renderer','Painters')
print('-depsc','-tiff','-r300','-painters',[figdir,filename1])

