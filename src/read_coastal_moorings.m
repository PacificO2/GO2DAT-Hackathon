function [lon,lat] = read_WOD
% This function reads the locations of all the coastal moorings off Chile.
% At the moment only the location.
%
% This script is needed by PacificO2_main.m
%
% History:
% Created the 3.May.2022 by the PacificO2 team for the GO2Dat Hackaton
% 2022.

% clear all; close all; clc

%% === Read directories and file names ===

% File containing the mooring information not yet published
moordir = '/Users/ana/Data/hackathon_go2dat_2022/Data/Raw/Coastal_moorings/';
fname = [moordir,'Coastal_moorings_location.xlsx'];

[NUM,TXT,RAW] = xlsread(fname);

lat = NUM(:,1);
lon = NUM(:,2);

end