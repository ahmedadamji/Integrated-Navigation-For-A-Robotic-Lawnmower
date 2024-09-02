%Solution

%% This script was written by Ikedinaekpere Kennedy Dike, Ahmed Adamjee,
%% and Abdulbaasit Sanusi

% This script calls all scripts that are required to fusing data from all 
% sensors to obtain an optimal position, velocity and heading solution.

% The following commands are used to clear the Command Window and
% Workspace.
clear;
close all;
clc;


%% Calling this function generates horizontal position and horizontal velocity 
% using GNSS only using an Extended Kalman Filter
GNSS_computation = GNSS_EKF_function;

%% figures showing the latitude, longitude and velocity of the lawnmower 
% using GNSS results with respect to time
gnss_latitude = GNSS_computation(:,2);
gnss_longitude = GNSS_computation(:,3);
figure
geoplot(gnss_latitude,gnss_longitude)
title('A plot of the position of the lawnmower in the map')
%clearing variables after plot
clearvars("gnss_longitude","gnss_latitude")


%% Calling this function generates horizontal position, horizontal velocity 
% and heading solution solution using all the other sensors apart from GNSS
% using an Extended Kalman Filter, however it uses the GNSS computation to
% initalise the Extended Kalman Filter
DR_computation = Multisensor_EKF_function(GNSS_computation);

%% figures showing the latitude, longitude and velocity of the lawnmower 
% using multisensor results with respect to time
DR_latitude = DR_computation(:,2);
DR_longitude = DR_computation(:,3);
figure
geoplot(DR_latitude,DR_longitude)
title('A plot of the position of the lawnmower in the map')
%clearing variables after plot
clearvars("DR_latitude","DR_longitude")


%% Calling this function generates horizontal position, horizontal velocity 
% and heading solution solution by combining the data from GNSS and the
% other sensors
Integration_result = Integration_function(GNSS_computation,DR_computation);

%% figures showing the latitude, longitude and velocity of the lawnmower 
% using GNSS as well as multisensor results with respect to time
Integration_latitude = Integration_result(:,2);
Integration_longitude = Integration_result(:,3);
figure
geoplot(Integration_latitude,Integration_longitude)
title('A plot of the position of the lawnmower in the map')
%clearing variables after plot
clearvars("Integration_latitude","Integration_longitude")
