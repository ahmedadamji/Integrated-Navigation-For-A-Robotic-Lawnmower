%Multisensor_EKF_function
function DR_computation = Multisensor_EKF_function(GNSS_computation)

%% This function was created on 08/02/2022 by Ikedinaekpere Kennedy Dike, Ahmed Adamjee,
%% and Abdulbaasit Sanusi

%% The aim of this function is to applying Kalman filtering to integrating
% data from different sensors
%
% Inputs:
%   GNSS_computation    continuous GNSS navigation solution
%     columns 1            Time (s)
%     columns 2            GNSS estimated user latitude (degree)
%     columns 3            GNSS estimated user longitude (degree)
%     columns 4            GNSS estimated user height (m)
%     columns 5-6          GNSS estimated user velocity (m/s)
%     columns 7            GNSS estimated user heading (degree)
%
% Outputs:
%   DR_computation      continuous Multisensor navigation solution
%     columns 1            Time (s)
%     columns 2            Multisensor estimated user latitude (degree)
%     columns 3            Multisensor estimated user longitude (degree)
%     columns 4-5          Multisensor estimated user damped velocity (m/s)
%     columns 6            Multisensor estimated user heading (degree)
%
%% Begins

%% Calling this script defines all the constants required for this Application.
Define_Constants_Master;

%% Computing the heading using the gyroscope's measured rate
gyro_heading(1) = gyro_ang_rate(1) * (time(2)-time(1));
for i = 2:length(time)
    gyro_heading(i) = gyro_heading(i-1) + (gyro_ang_rate(i) * (time(i)-time(i-1)));
end
correct_gyro_heading_list(1) = gyro_heading(1);

%% Initialising state estimates, Error Covariance Matrix for the gyroscope
% readings, to update using an extended Kalman Filter
gyro_filter_states = zeros(2,1);
P = [gryo_noise_std 0; 0 gyro_bias];
% tau is the propagation delay 
tau = 0.5;

%% Using the transition matrix to propagate the state estimates for the gyroscope readings %%
% Creating a loop to do this for all time epochs
for i = 2 : length(time)
    % Transition matrix (Φ k-1)
    transition_matrix = [1 tau;0 1];

    % Computing the system noise covariance matrix 
    Q_k_1 = [(gyro_measurement_errors*tau)+((1/3)*gyro_measurement_errors*tau^3), ((1/2)*gyro_measurement_errors*tau^2);...
        ((1/2)*gyro_measurement_errors*tau^2),gyro_measurement_errors*tau ];

    % Using the transition matrix to propagate the state estimates for 
    % the gyroscope readings for the current epoch
    gyro_filter_states = transition_matrix*gyro_filter_states;

    % propagating the error covariance matrix 
    P = (transition_matrix*P*transition_matrix') + Q_k_1;

    % Computing the measurement matrix (H)
    H = [-1 0];

    % Computing the measurement noise covariance matrix (R)
    R_k = magnetic_heading_noise_std;

    % Computing the Kalman gain matrix (K)
    K_k = (P* H') /((H*P*H')+R_k);

    % Formulating the measurement innovation vector (delta_z)
    delta_z = (heading(i)-gyro_heading(i)) - (H*gyro_filter_states);

    % Updating the state estimates %
    gyro_filter_states = gyro_filter_states + (K_k*delta_z);

    % Updating the error covariance matrix %%
    P = (eye(2) - (K_k*H))*P;

    % Updating correct gyroscope heading values
    correct_gyro_heading_list(i) = gyro_heading(i) - gyro_filter_states(1);

end

%% Ensuring the heading readings are between pi and -pi
% correct_gyro_heading_list(correct_gyro_heading_list>pi) = correct_gyro_heading_list(correct_gyro_heading_list>pi) - (2*pi);
% correct_gyro_heading_list(correct_gyro_heading_list<(-pi)) = correct_gyro_heading_list(correct_gyro_heading_list<(-pi)) + (2*pi);


%% Computing forward speed from sensors for each wheel by computing their average
% Please note that this is a simple estimate
forward_speed = mean(forward_speed,2);%Ask TA what to do here. Mean it or not

% The initial geodetic latitude and longitude (at time 0) are 
% used from the respective corresponding GNSS solution
L_k =  GNSS_computation(1,2) * deg_to_rad;
lamda_k = GNSS_computation(1,3) * deg_to_rad;
% The height solution at all epochs are used from the respective 
% corresponding GNSS solution
height_list = GNSS_computation(:,4);
% Initialising heading solution from previous epoch
correct_gyro_heading_k_1 = correct_gyro_heading_list(1) ;
% Initialising heading solution from current epoch
correct_gyro_heading_k = correct_gyro_heading_list(2);

% instantaneous velocity at time zero is given by the speed and heading
% measurements at that time:
v_damped = [(forward_speed(1)*(cos(correct_gyro_heading_k_1)));(forward_speed(1)*(sin(correct_gyro_heading_k_1)))];
% Initialising a list that stores all the damped instantaneous DR velocity
% at each epoch:
v_damped_list = [v_damped];

% Storing the longitude and longitude in degrees to array of results
longitude_list = [lamda_k * rad_to_deg];
latitude_list = [L_k * rad_to_deg];

%% Use the speed and heading measurements to compute a position solution for
% the rest of the car s trajectory. %%

for k=2:length(time)
    correct_gyro_heading_k_1 = correct_gyro_heading_list(k-1);
    correct_gyro_heading_k = correct_gyro_heading_list(k);
    % average velocity between epochs k−1 and k is given by
    v = 0.5* [(cos(correct_gyro_heading_k) + cos(correct_gyro_heading_k_1));(sin(correct_gyro_heading_k) + sin(correct_gyro_heading_k_1))]*forward_speed(k);

    %% RN (meridian radius of curvature) and RE (transverse radius of 
    % curvature) are computed from the latitude using Radii_of_curvature().
    [R_N,R_E]= Radii_of_curvature(L_k);
    
    %% h is the geodetic height, computed using the GNSS solution
    h = height_list(k);
    % ^^ I HAVE CHANGED THIS FROM 1 to K, please confirm!!!!!

    %% The latitude, and longitude, at epoch k are computed from using values
    % at epoch k − 1 using.%%
    L_k = L_k + (v(1)*(time(k) - time(k-1)))/(R_N+h);
    lamda_k = lamda_k + (v(2)*(time(k) - time(k-1)))/((R_E+h)*cos(L_k));


    %% Computing the damped instantaneous DR velocity at each epoch %%
    v_damped(1) = (1.7*v(1)) - (0.7*v_damped(1));
    v_damped(2) = (1.7*v(2)) - (0.7*v_damped(2));

    %% Saving the results to a list %%
    v_damped_list = [v_damped_list v_damped];
    longitude_list = [longitude_list; lamda_k * rad_to_deg];
    latitude_list = [latitude_list; L_k * rad_to_deg];

end
%transposing v_damed_list to correct shape to store in results array
v_damped_list = v_damped_list';

%% Storing all results from the Multisensor Extended Kalman Filter
% solution in a matrix called DR_computation
DR_computation = [time latitude_list longitude_list v_damped_list correct_gyro_heading_list'*rad_to_deg];
end