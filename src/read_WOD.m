function [lon,lat] = read_WOD
% This function reads the WOD database, previously downloaded from:
% https://www.ncei.noaa.gov/access/world-ocean-database-select/dbsearch.html
% The data read here was already subsampled in the WOD viewer and contains
% data from the region -90 to -65 E; and -40 to 10 N.
%
% This script is read by PacificO2_main.m
%
% History:
% Created the 3.May.2022 by the PacificO2 team for the GO2Dat Hackaton
% 2022.

% clear all; close all; clc
%% === Read directories and file names ===

% All the files from the Mirai cruise were downloaded to this directory:
woddir = '/Users/ana/Data/hackathon_go2dat_2022/Data/Raw/WOD_PacificO2/wod_ncfiles/';

% This reads the names of all the files that end in .csv, which contain
% the data that intersts us.
fnames = dir([woddir,'*.nc']);

% === Read the data from each file ===

% Create the variables to be read in the loop
lat = []; lon = []; t_utc = []; pressure = [];
oxy = []; oxy_flag = []; temp = []; temp_flag = [];
salt = []; salt_flag = []; wod_unique_cast = [];

% This loop that opens each individual csv file, finds the temperature,
% salinity, oxygen, lat and lon and depth (in this case, intake depth) for
% each file
for i = 1:length(fnames)
    fname = [woddir fnames(i).name]; % It has troubles reading the EXPOCODE, so start reading from the second column
    
    %     t_utc = cat(1,t_utc,ncread(fname,'GMT_time'));
    lat = cat(1,lat,ncread(fname,'lat'));
    lon = cat(1,lon,ncread(fname,'lon'));
    %     pressure = cat(1,pressure,ncread(fname,'z')); % In reality this is depth (in m); This would need to be converted to Pressure. or not.
    %
    %     wod_unique_cast = cat(1,wod_unique_cast,ncread(fname,'Cast_Tow_number'));
    %
    %     oxy = cat(1,oxy,ncread(fname,'Oxygen')); % In µmol/kg
    %     oxy_flag = cat(1,oxy_flag,ncread(fname,'Oxygen_WODflag')); %
    %
    %     temp = cat(1,temp,ncread(fname,'Temperature')); % "degree_C"
    %     temp_flag = cat(1,temp_flag,ncread(fname,'Temperature_WODflag'));
    %
    %     salt = cat(1,salt,ncread(fname,'Salinity')); % "degree_C"
    %     salt_flag = cat(1,salt_flag,ncread(fname,'Salinity_WODflag'));
    
end

% === Questions:
% - How to link the coordinates with each profile?
end