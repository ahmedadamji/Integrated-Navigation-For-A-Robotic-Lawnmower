%leastsquares
function [x_est,P_matrix] = leastsquares
%% This function created 26/01/2022 by Ahmed Adamjee, Kenndy Dike and
%% Abdulbaasit Sanusi
%
%% The aim of this function script is to initialise
% the GNSS EKF state estimates and error covariance matrix using
% Least squares approach for Coursework 1.

%
% Outputs:
%   x_est                 Kalman filter initial estimates:
%     Rows 1-3            estimated ECEF user position (m)
%     Rows 4-6            estimated ECEF user velocity (m/s)
%     Row 7               estimated receiver clock offset (m)
%     Row 8               estimated receiver clock drift (m/s)

%   P_matrix              state estimation error covariance matrix

%% Begins

%% Obtain constant parameters and files needed
Define_Constants_Master;

ranges = Pseudoranges(2,2:end);
range_rates = Pseudorangerates(2,2:end);

%% Initial position and velocity estimates
r_ea = [0;0;0];
v_ea = [0;0;0];


%% Compute the Cartesian ECF position and velocity of the satelittes at time 0
r_ej = cell(length(satellite),1);
sat_v_es_e = cell(length(satellite),1);
for i=1:length(satellite)
    [r_ej{i,1},sat_v_es_e{i,1}] = Satellite_position_and_velocity(0,satellite(i));
    % positions and velocity for each satelitte is stored in a cell
end

%% Initialise clock offset and drift errors
pc =0;
drift = 0;

%% Form initial state estimate
x = [r_ea; pc];
xdot = [v_ea; drift];

%% Least sqaure Algorithm begins. Iterate till convergence
for k = 1:4

    r_ea = x(1:3); pc = x(4); 
    v_ea = xdot(1:3); drift = xdot(4);
    
    %% Predict ranges and range rates from the approx. user positiona and velcoity to each satelitte 
    sagnac_effect = cell(length(satellite),1);
    r_aj = [];
    r_aj_dot = [];
    u_aj = cell(length(satellite),1);
    r=0;
    rdot = 0;

   for i=1:length(satellite)
        % The recursion is resolved by initially computing the range with
        % the Sagnac effect compensation matrix set to the identity matrix
        C_e{i,1} = eye(3,3);
        % computimng the corresponf=ding value for the range and range
        % rate
        r = sqrt(transpose ( (C_e{i,1}*transpose(r_ej{i,1})) -r_ea)* ((C_e{i,1}*transpose(r_ej{i,1}))-r_ea));
        first_term = C_e{i,1} * (sat_v_es_e{i,1}'+(Omega_ie*r_ej{i,1}'));
        second_term = v_ea +  (Omega_ie*r_ea);
        r_dot =  (first_term - second_term);

        % using this range to compute the Sagnac effect compensation matrix
        C_e{i,1} = [1 (omega_ie*r/c) 0;(-omega_ie*r/c) 1 0;0 0 1]; % Sagnac effect compensation matrix.
        % Recomputing the range and range rates.
        r = sqrt(transpose ( (C_e{i,1}*transpose(r_ej{i,1})) -r_ea)* ((C_e{i,1}*transpose(r_ej{i,1}))-r_ea));
        first_term = C_e{i,1} * (sat_v_es_e{i,1}'+(Omega_ie*r_ej{i,1}'));
        second_term = v_ea +  (Omega_ie*r_ea);
        r_dot =  (first_term - second_term);

        % Computing the line-of-sight unit vector from the approximate user
        % position to each satellite
        u_aj{i,1} = (((C_e{i,1}*transpose(r_ej{i,1}))-r_ea))/r;
        % Predicting the range rates from the approximate user position to each satellite
        r_dot = u_aj{i,1}' * r_dot;

        % Appending the range and range rates to an array to store result
        % for each satellite
        r_aj = [r_aj r];
        r_aj_dot = [r_aj_dot r_dot];

    end
    
    %% Compute Measurement Matrix
    H = [];
    for i = 1:length(satellite)
        H = [H;-u_aj{i,1}' 1];
    end
    
    %% Compute Measurement Innovation Vector
    delta_z = [];
    delta_z_dot = [];
    for i = 1:length(satellite)
        delta_z = [delta_z (ranges(i)-r_aj(i)-pc)];
        delta_z_dot = [delta_z_dot (range_rates(i)-r_aj_dot(i)-drift)];

    end
    delta_z = transpose(delta_z);
    delta_z_dot = transpose(delta_z_dot);
    
    %% Compute the position and receiver clock offset
    x = x + (inv(transpose(H)*H) * transpose(H)* delta_z);
    
    %% Compute the velcoity and receiver clock drift
    xdot = xdot + (inv(transpose(H)*H) * transpose(H)* delta_z_dot);

end

%% Convert Cartesian ECEF positiona and vecolity solution to longitude, latitude, height and velocity
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x(1:3),xdot(1:3));

L_b = rad2deg(L_b); % Convert to degrees
lambda_b = rad2deg(lambda_b); % Convert to degrees

%% Initialise state estimates to use fo Kalman Filter
x_est = [x(1:3);xdot(1:3);x(4);xdot(4)];

%% Initialise error covariance matrix
P_matrix =  zeros(8);
P_matrix(1,1) = (noise_std_on_pseudo_range^2);
P_matrix(2,2) = (noise_std_on_pseudo_range^2);
P_matrix(3,3) = (noise_std_on_pseudo_range^2);
P_matrix(4,4) = (noise_std_on_pseudo_rangerates^2);
P_matrix(5,5) = (noise_std_on_pseudo_rangerates^2);
P_matrix(6,6) = (noise_std_on_pseudo_rangerates^2);
P_matrix(7,7) = (noise_std_on_pseudo_range^2);
P_matrix(8,8) = (noise_std_on_pseudo_rangerates^2);


%% Ends










