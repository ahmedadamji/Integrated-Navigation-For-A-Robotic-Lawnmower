%GNSS_EKF_function
function GNSS_computation = GNSS_EKF_function

%% This function was created on 08/02/2022 by Ikedinaekpere Kennedy Dike, Ahmed Adamjee,
%% and Abdulbaasit Sanusi

%% The aim of this function is to apply Kalman filtering to compute a
% continuous GNSS navigation solution using an initial estimate from
% least-squares estimation
%
% Outputs:
%   GNSS_computation    continuous GNSS navigation solution
%     columns 1            Time (s)
%     columns 2            GNSS estimated user latitude (degree)
%     columns 3            GNSS estimated user longitude (degree)
%     columns 4            GNSS estimated user height (m)
%     columns 5-6          GNSS estimated user velocity (m/s)
%     columns 7            GNSS estimated user heading (degree)
%
%% Begins

%% Calling this script defines all the constants required for this Application.
Define_Constants_Master;

%% Initalising arrays and variables used in this script

% Boolean variable that is set to true if any position solution found is an
% outlier
outlier_found = false;
% integer to store how many outliers were found in current epoch
num_outliers = 0;
% integer to store the number of stellites that did not have outliers as
% this is needed to determine if we need to run outlier detection again
satellites_without_outliers = 0;

% List to store velocity solution from each epoch
velocity_list = [];
% List to store height solution from each epoch
height_list = [];
% List to store latitude solution from each epoch
latitude_list = [];
% List to store longitude solution from each epoch
longitude_list = [];
% List to store heading solution from each epoch
heading_list = [];

%% initialising state estimates and the error covariance matrix %%

% Calling this script initialises our state estimates and the error
% covariance matrix
[x_est, P_matrix] = leastsquares;
x_est;
P_matrix;

% Decomposing the state vector estimate into position,
% velocity, clock offset and clock drift
r_ea = x_est(1:3);
v_ea = x_est(4:6);
offset = x_est(7);
drift = x_est(8);

% Renaming the error covariance matrix to P
P = P_matrix;

%% Computing the transition matrix %%

% Setting the propagation delay to 1
propagation_delay = 0.5; %tau_s %%NOT SURE AS IT IS NOT GIVEN

% Identity matrix of size 3x3
I = eye(3,3);
% Matrix of zeros of size 3x3
O = zeros(3,3);

% Transition matrix (Φ k-1)
transition_matrix = [I (propagation_delay*I) zeros(3,1), zeros(3,1);...
    O I zeros(3,1), zeros(3,1);...
    zeros(1,3) zeros(1,3) 1 propagation_delay;...
    zeros(1,3) zeros(1,3) 0 1];


%% Computing the system noise covariance matrix %%

Q_k_1 = zeros(8,8); % initialising the system noise covariance matrix
% CHECK THIS VALUE -->
S_a = 0.05; % the accele0* ration power spectral density (PSD), is 5 m2s−3. %NOT SURE AS IT IS NOT GIVEN
S_c_phi = clock_phase_PSD; % the clock phase PSD, is 0.01 m2s−3
S_c_f = clock_freq_PSD; % the clock frequency PSD, is 0.04 m2s−1
Q_k_1(1:3,1:3) = (1/3)*S_a*power(propagation_delay,3)*I;
Q_k_1(1:3,4:6) = (1/2)*S_a*power(propagation_delay,2)*I;
Q_k_1(1:3,7) = zeros(3,1);
Q_k_1(1:3,8) = zeros(3,1);
Q_k_1(4:6,1:3) = (1/2)*S_a*power(propagation_delay,2)*I;
Q_k_1(4:6,4:6) = S_a*propagation_delay*I;
Q_k_1(4:6,7) = zeros(3,1);
Q_k_1(4:6,8) = zeros(3,1);
Q_k_1(7,1:3) = zeros(1,3);
Q_k_1(7,4:6) = zeros(1,3);
Q_k_1(7,7) = (S_c_phi*propagation_delay) + ((1/3)*S_c_f*power(propagation_delay,3));
Q_k_1(7,8) = ((1/2)*S_c_f*power(propagation_delay,2));
Q_k_1(8,1:3) = zeros(1,3);
Q_k_1(8,4:6) = zeros(1,3);
Q_k_1(8,7) = ((1/2)*S_c_f*power(propagation_delay,2));
Q_k_1(8,8) = (S_c_f*propagation_delay);


%% Using the transition matrix to propagate the state estimates %%
% Creating a loop to do this for all time epochs
for k = 2: length(time)+1

    % Array to stores satellite range data for current epoch
    ranges = Pseudoranges(k,2:end);
    % Array to stores satellite range rates data for current epoch
    range_rates = Pseudorangerates(k,2:end);

    % Using the transition matrix to propagate the state estimates
    x_est = transition_matrix*x_est;

    % Decomposing the state vector estimate into position,
    % velocity, clock offset and clock drift
    r_ea = x_est(1:3);
    v_ea = x_est(4:6);
    offset = x_est(7);
    drift = x_est(8);

    %% propagating the error covariance matrix %%
    P = (transition_matrix*P*transpose(transition_matrix))+Q_k_1;

    %% Predicting the ranges from the approximate user position to each satellite %%

    % Initialising an array to store the ranges from the approximate user position to each satellite
    r_aj = [];
    % Initialising cells to store the Cartesian ECEF position and velocity of satellite
    r_ej = cell(length(satellite),1);
    sat_v_es_e = cell(length(satellite),1);
    for i=1:length(satellite)
        % Computing the Cartesian ECEF position and velocity for each satellite
        % using the Matlab function Satellite_position_and_velocity.m
        [r_ej{i,1},sat_v_es_e{i,1}] = Satellite_position_and_velocity(Pseudoranges(k,1),satellite(i));
    end
    % Initialising cells to store the Sagnac effect compensation matrix for
    % each satellite
    C_e = cell(length(satellite),1);
    % Initialising cells to store the line-of-sight unit vector (u_aj) from the
    % approximate user position to each satellite:
    u_aj = cell(length(satellite),1);
    % r stores the ranges from the approximate user position for current
    % satellite
    r=0;
    % r_dot stores the range rates from the approximate user position for
    % current satellite
    r_dot = 0;
    % r_dot stores the range rates from the approximate user position to each satellite
    r_aj_dot = [];

    %Itterating through each satellite to compute ranges and range rates
    % to each satellite
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

    %% Computing the measurement matrix (H) %%
    H = [];
    for i = 1:length(satellite)
        H = [H; -(u_aj{i,1})' 0 0 0 1 0];
    end
    for i = 1:length(satellite)
        H = [H;0 0 0 -(u_aj{i,1})' 0 1];
    end


    %% Computing the measurement noise covariance matrix (R) %%
    % sigma_p --> error standard deviation on pseudorange measurements
    sigma_p = 10^2;
    % sigma_r --> error standard deviation on pseudo-range rate measurements
    sigma_r = 0.05^2;
    n = length(H)/2;
    R = [sigma_p*(eye(n)), zeros(n); zeros(n), sigma_r*(eye(n))];


    %% Computing the Kalman gain matrix (K) %%
    K = P*transpose(H)/(( H*P*transpose(H)) + R);


    %% Formulating the measurement innovation vector (delta_z) %%
    delta_z = [];
    for i = 1:length(satellite)
        % ranges(i) is the measured pseudo-range from satellite i to the 
        % user antenna
        % offset is the propagated receiver clock offset estimate
        delta_z = [delta_z ;(ranges(i)-r_aj(i)-offset)];
    end
    for i = 1:length(satellite)
        % range_rates(i) is the measured pseudo-range rate from satellite i
        % to the user antenna
        % drift s the propagated receiver clock drift estimate
        delta_z = [delta_z ;(range_rates(i)-r_aj_dot(i)-drift)];
    end

    %% --------------------------------------------------------------------
    %% A Residual-based outlier detection at each epoch.


    % Initialising an identity matrix, used to Compute outliers
    I_m = eye(16,16);

    % Computing the residuals vector (v)
    v = (H* inv(H' *  H) * H' - I_m) * delta_z;

    % Computing the residuals covariance matrix
    % Im is the mxm identity matrix, where m is the number of measurements
    % noise_std_on_pseudo_range is the measurement error standard deviation
    C_v = (I_m - H* inv(H' *  H) * H')* noise_std_on_pseudo_range^2;

%     Computing the normalised residuals and comparing
%     them for each satellite with a threshold
    for i = 1:length(satellite)
        if (k>2) %To start finding solution from epoch 2 onwards
            % Measurement from satellite i is an outlier when the following
            % condition is met:
            if norm(v(i)) > (sqrt(C_v(i,i))*outlier_threshold)
                % Setting the result to the boolean variable outlier_found
                outlier_found = true;


                %disp(time(k-1));
                %disp(satellite(i));
                
                % Now excluding this measurement from our calculations 
                % by deleting the corresponding data from our calculations
                % for both range and range rates
                delta_z(i) = [];
                delta_z(length(satellite)+i) =[];
                H(i,:) = [];
                H(length(satellite)+i,:)=[];
                K(:,i) = [];
                K(:,length(satellite)+i) = [];
                %increasing the total number of outliers found in
                %num_outliers
                num_outliers = num_outliers+1;
                break;
            end
        end
    end

    % If outliers were detected, we recalculate our position at that epoch 
    % without the measurement that had the largest residual, retaining any
    % outlying measurements with smaller residuals
    while(outlier_found)

        % Initialising an identity matrix, used to Compute outliers
        I_m = eye(16-(num_outliers*2),16-(num_outliers*2));
        
        % Computing the residuals vector (v)
        v = (H* inv(H' *  H) * H' - I_m) * delta_z;

        % Computing the residuals covariance matrix
        % Im is the mxm identity matrix, where m is the number of measurements
        % noise_std_on_pseudo_range is the measurement error standard deviation
        C_v = (I_m - H* inv(H' *  H) * H')* noise_std_on_pseudo_range^2;

        % Computing the normalised residuals and comparing
        % them for each satellite with a threshold
        for i = 1:length(satellite)-num_outliers % -num_outliers as we have already eliminated readings from some satellites
            if (k>2) %To start finding solution from epoch 2 onwards
                % Measurement from satellite i is an outlier when the following
                % condition is met:
                if norm(v(i)) > (sqrt(C_v(i,i))*outlier_threshold)
                    % disp(time(k-1));
                    % disp(satellite(i));
                    % The following message is displayed if an outlier was
                    % yet detected even after eliminating the largest
                    % outlier, in this case we must create another loop to
                    % recompute if we still have more outliers
                    disp("An outlier was found, please" + ...
                        "check if still remains!");

                    %disp(num_outliers);

                    % Now excluding this measurement from our calculations 
                    % by deleting the corresponding data from our calculations
                    % for both range and range rates
                    delta_z(i) = [];
                    delta_z(length(satellite)-num_outliers+i) =[];
                    H(i,:) = [];
                    H(length(satellite)-num_outliers+i,:)=[];
                    K(:,i) = [];
                    K(:,length(satellite)-num_outliers+i) = [];

                    %increasing the total number of outliers found in
                    %num_outliers
                    num_outliers = num_outliers+1;
                    break;
                else
                    satellites_without_outliers = satellites_without_outliers+1;
                end
            end
        end
        if (satellites_without_outliers==(length(satellite)-num_outliers))
            break;
        end
    end
    % setting outlier_found to false as the outliers have now been
    % eliminated
    outlier_found = false;
    num_outliers = 0;
    satellites_without_outliers=0;

    %% --------------------------------------------------------------------

    %% Updating the state estimates %%
    x_est = x_est + (K*delta_z);

    %% Updating the error covariance matrix %%
    P = (eye(8)-(K*H))*P;

    %% Converting this Cartesian ECEF position solution to latitude, 
    % longitude and height using the Matlab function pv_ECEF_to_NED.m 
    % This also converts the velocity solution from ECEF resolving axes 
    % to north, east and down %%
    [latitude,longitude,height,velocity] = pv_ECEF_to_NED(x_est(1:3),x_est(4:6));
    
    % the Matlab function outputs latitude and longitude in radians, so
    % converting to degrees:
    latitude = rad2deg(latitude);
    longitude = rad2deg(longitude);

    % We cannot find the heading using GNSS as it is only obtainable from
    % Magnetic compass, Trajectory (from positioning fixing) and
    % Gyroscope (After initialisation).
    heading = 0;

    %% Storing solutions found in each epoch to a list

    latitude_list = [latitude_list, latitude];

    longitude_list = [longitude_list, longitude];

    height_list = [height_list, height];

    velocity_list = [velocity_list, velocity];

    heading_list = [heading_list, heading];

end

%% Storing all results from the GNSS Extended Kalman Filter
% solution in a matrix called GNSS_computation
GNSS_computation =[time,latitude_list',longitude_list',height_list',velocity_list(1:2,:)', heading_list'];


% plot(time,longitude_list)
% title('A plot of the longitude of the lawnmower with respect to time')

% figure
% 
% plot(time,velocity_list(1,:))
% title('A plot of the velcoity in the east direction of the lawnmower with respect to time')
% figure
% 
% plot(time,velocity_list(2,:))
% title('A plot of the velcoity in the north direction of the lawnmower with respect to time')
end