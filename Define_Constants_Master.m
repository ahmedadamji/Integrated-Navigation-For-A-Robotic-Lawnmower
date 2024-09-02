%Define_Constants_Master

%% This script was written by Ikedinaekpere Kennedy Dike, Ahmed Adamjee,
%% Abdulbaasit Sanusi

%% This script defines all the constants required for this application. It also
%% makes use of the "Define_Constants.m" File which contains some constants
%% required
Define_Constants;

%% Reading Pseudorange and Pseudorangerates data
Pseudoranges = importdata("Pseudo_ranges.csv");
Pseudorangerates = importdata("Pseudo_range_rates.csv");
%Array storing all time values
time = (Pseudoranges(2:end,1));
%Array storing all satellite numbers
satellite = Pseudoranges(1,2:end);

%% Reading Dead reckoning data
Speed_Heading = importdata("Dead_reckoning.csv");
%Array storing all time values
time = Speed_Heading(:,1);
%Array storing the wheel-speed measurements from sensors 1 to 4 in 
% metres per second
forward_speed = Speed_Heading(:,2:5); %Ask in lab how many wheels and if we should take average speed at centre
%Array storing the measurements in degrees from the magnetic compass
gyro_ang_rate = Speed_Heading(:,6);
%Array storing the heading in degrees.
heading = Speed_Heading(:,7)*deg_to_rad;

%% Sensor Specifications

%% GNSS
signal_space_error_std = 1;
residual_ionosphere_error_std = 2;
residual_troposphere_error_std = 0.2;
code_tracking_multipath_error_std = 2;
range_rate_tracking_multipath_error_std = 0.02;
noise_std_on_pseudo_range = 10; %used in outlier detection however in WS1 it was called measurement error std confirm if it's the same
noise_std_on_pseudo_rangerates = 0.05;
outlier_threshold = 4;  %Not given BUT 6 was given in the WS1
reciver_clock_offset_std = 100000;
reciver_clock_drift_std = 200;
clock_phase_PSD = 0.01;
clock_freq_PSD = 0.04;

%% Location of wheel speed sensors
wheel_sensor1 = [0.3 -0.2]; %[Forwards right]
wheel_sensor2 = [0.3 0.2]; %[Forwards right]
wheel_sensor3 = [-0.2 -0.2]; %[Forwards right]
wheel_sensor4 = [-0.2 0.2]; %[Forwards right]

% The Following are assumed to account for the effects of varying tyre 
% radii and wheel slip
wheel_scale_factor_error_std = 0.03; %3 percent
wheel_noise_std = 0.05;
wheel_speed_sensor_outputs = 0.02;

%% Gyroscope Error Characteriscts
gyro_bias = 1 * deg_to_rad;
gyro_scale_factor_error = 0.01; % 1 percent
gryo_cross_coupling_error = 0.1/100;
gryo_noise_std = 10^(-4);
gryo_quant_level = 2*10^(-4);

power_spectral_desnity = 0.01; %wheel speed errors
gyro_measurement_errors = 3*10^(-6); %Does not include bias
velocity_error_variance_growth = 0.01;
heading_error_variance_growth = 3*10^(-6);

magnetic_heading_noise_std = 4 * deg_to_rad;% Does not consider absence of localised magnetic anomalies.
