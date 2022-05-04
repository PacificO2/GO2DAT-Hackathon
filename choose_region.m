% GO2DAT-Hackathon

% This script reads the MIRAI database, previously downloaded from: 
% https://drive.google.com/file/d/1lbh7zfDA_jznyQ6yJ19gDgMhrBLttZgs/view
% Once it reads the data and files, it selects a box with the desired
% coordinates. 
% Thi script is read by PacificO2_main.m
% 
% History: 
% Created the 3.May.2022 by the PacificO2 team for the GO2Dat Hackaton
% 2022. 

clear all; close all; clc

%% === Read directories and file names ===

% All the files from the Mirai cruise were downloaded to this directory: 
% miraidir = '/Users/ana/Data/hackathon_go2dat_2022/Data/Raw/mirai_tsg_dataset_20220418/';
miraidir = "C:\Users\user\Documents\MATLAB\mirai_tsg_dataset_20220418";

% This reads the names of all the files that end in .csv, which contain
% the data that intersts us. 
fnames = dir(['49N','*.csv']);

% === Read the data from each file ===

% Create the variables to be read in the loop
lat = []; lon = []; t_utc = []; pressure = [];
oxy = []; oxy_flag = []; temp = []; temp_flag = []; 
salt = []; salt_flag = []; 

% This loop that opens each individual csv file, finds the temperature,
% salinity, oxygen, lat and lon and depth (in this case, intake depth) for
% each file
for i = 1:length(fnames)
    fname = [miraidir fnames(i).name]; % It has troubles reading the EXPOCODE, so start reading from the second column
    tmp = csvread(fnames(i).name,3,1);
    
    % Find the variables, it would take too much time to make matlab detect
    % the correct column automatically, we had to ASSUME that the columns
    % for each variable do not change. 
    %t_utc = cat(1,t_utc(:,2)); % Time UTC
    lat = cat(1,lat,tmp(:,3)); % Concatenate the latitude
    lon = cat(1,lon,tmp(:,4)); % Concatenate the longitude
    pressure = cat(1,pressure,tmp(:,6)); % Concatenate the pressure 
    
    oxy = cat(1,oxy,tmp(:,13)); % I selected the variable named TSG-OXYGEN (in µmol/kg)
    oxy_flag = cat(1,oxy_flag,tmp(:,14)); % I selected the variable named TSG-OXYGEN (in µmol/kg)
    
    temp = cat(1,temp,tmp(:,7)); % Variable named INTAKE-TEMPERATURE (ITS-90)
    temp_flag = cat(1,temp_flag,tmp(:,8)); 
    
    salt = cat(1,salt,tmp(:,9)); % Variable named TSG-SALINITY (PSS-78)
    salt_flag = cat(1,salt_flag,tmp(:,10));
end
    
% Questions for MIRAI: 
% - what does it mean TSG-OXYGEN? 
% - Is the data contained in another database? 
   
%% === Select a box with the desired coordinates ===

% edit and choose desired coordinates 
% note that you would have error if you don't choose desired coordinates and
%data
lower_lon = []; 
higher_lon = [];
lower_lat = [];
higher_lat = []; 
data = [] %choose desired data e.g oxy,salt,temp.

%rearrage into lon,lat and the desired data
lon_lat_data = cat(2,lon,lat,data);

% index to subset to desired region
ind = find(lat >= lower_lat & lat<= higher_lat & lon >= lower_lon & lat<= higher_lon);

new_data = data(ind);
new_lon = lon(ind);
new_lat = lat(ind);

% data for desired region
for i = 1:length(new_lon)
    desired_data(i,1) = new_lon(i);
    desired_data(i,2) = new_lat(i);
    desired_data(i,3) = new_data(i);
end

%output is the variable 'desired_data' containing lon, lat and the data