% This script compiles O2 data from a box off Chile, Peru and Ecuador with
% delimieted by: -90 to -65 E; and 10 to 40 S. 
% The databases compiled are: 

% 1) Ocean Data base --> read_WOD.m

% 2) Mirai Databse --> read_Mirai.m

% 3) Coastal moorings --> read_coastal_moorings

% This main script compiles the data into one file (for the moment just
% .mat) that can then be used with the data_explotaition "toolbox". 

% History:
% Created the 3.May.2022 by the PacificO2 team for the GO2Dat Hackaton
% 2022.

clear all; close all; clc

%% === User defined ===

plot_map = true; 

%% === Compile all the databases regardless of instrument ===
% This would be like reading a GO2DAT product level 4

db = []; % This variable saves the database where the original data is comming from. 
% === Read the open ocean databases ===
[lon1,lat1] = read_WOD;
str = 'WOD';
db = cat(1,db,repmat({str},[size(lon1),1])); clear str

[lon2,lat2] = read_Mirai;
str = 'Mirai';
db = cat(1,db,repmat({str},[size(lon2),1])); clear str

% === Read the coastal databases ===
[lon3,lat3] = read_coastal_moorings; 
str = 'Coast_mooring';
db = cat(1,db,repmat({str},[size(lon3),1])); clear str

%% === Concatenate the information ===

lon = cat(1,lon1,lon2,lon3); %clear lon1 lon2 lon3
lat = cat(1,lat1,lat2,lat3); %clear lat1 lat2 lat3

%% === Make a plot of the available data ===
if plot_map
    m_proj('mercator','longitudes',[-90 -65],'latitudes',[-40 10]);
    m_gshhs_h('save','gumby');
    
    
    f1 = figure;
    m_grid('tickdir','out')
    hold on
    h = m_scatter(lon1,lat1,'k.'); % Plot the WOD datapoints
    hold on
    h = m_scatter(lon2,lat2,'bo'); % Plot the Mirai datapoints
    hold on
    h = m_scatter(lon3,lat3,100,'ro','filled'); % Plot the Coastal moorings
    hold on
    m_usercoast('gumby','patch',[0.8 0.8 0.8],'edgecolor','none');
    
    % m_legend('WOD','Mirai','coastal moorings')
    
    % === Save the figure ===
    figdir = '/Users/ana/Data/hackathon_go2dat_2022/figures/';
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
    
    filename1 = ['compiled_wod_mirai_coastal'];
    % for eps
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300','-painters',[figdir,filename1])
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save to .mat-file:

d.lon = lon;
d.lat = lat;
d.db = db;

d.filedetails.creationdate = datestr(now);
n = mfilename('fullpath');
d.filedetails.createdwith = [n,'.m'];
d.filedetails.compileddby = 'Ana C. Franco (afranco@eoas.ubc.ca)';
d.description = 'Data from several databases indicated in the variable db';

outputdir = '/Users/ana/Data/hackathon_go2dat_2022/Data/Processed/';
outputfile = [outputdir,'PacificO2_wod_mirai_coastal.mat'];

save(outputfile,'d');

