# Integrated Navigation for a Robotic Lawnmower

## Overview
This project implements sensor fusion techniques to estimate the horizontal position, velocity, and heading of a robotic lawnmower using data from GNSS, wheel speed sensors, a magnetic compass, and a MEMS gyroscope.

## Folder Structure
- `*.m`: MATLAB scripts for sensor data processing and estimation.
  - `GNSS_EKF_function.m`: Extended Kalman Filter implementation for GNSS data.
  - `Multisensor_EKF_function.m`: Multi-sensor fusion using EKF.
  - `Integration_function.m`: Combines data from different sensors.
  - `Satellite_position_and_velocity.m`: Computes satellite positions and velocities.
- `*.csv`: Input and example output data files.
  - `Dead_reckoning.csv`: Wheel speed sensor data.
  - `Pseudo_ranges.csv`: GNSS pseudo-range measurements.
  - `Pseudo_range_rates.csv`: GNSS pseudo-range rate measurements.
  - `Example_Output_Profile.csv`: Example output format.
- `README.md`: Project overview and usage instructions.
- `*.asv`: Autosave versions of MATLAB scripts.

## Task
- Correct GNSS measurements using least-squares for receiver clock errors.
- Implement an Extended Kalman Filter (EKF) for sensor fusion.
- Estimate the best possible position, velocity, and heading for the lawnmower.
- Generate a CSV file with time, latitude, longitude, north velocity, east velocity, and heading.

## Usage
1. Place the required `.csv` files in the project directory.
2. Run the main MATLAB scripts to process the data and compute estimates:
   - Start with `GNSS_EKF_function.m` to process GNSS data.
   - Use `Multisensor_EKF_function.m` for multi-sensor fusion.
3. The final motion profile will be generated in a CSV format.

## Requirements
- MATLAB or compatible environment.

## Notes
- The `Define_Constants.m` script is essential for setting up required constants.
- Utility scripts like `ECEF_to_NED.m` and `NED_to_ECEF.m` are provided for coordinate transformations.

## Acknowledgements
This project is based on the COMP0130 Robot Vision and Navigation coursework at UCL.
