% Clean up the workspace
clear all; close all; clc;

% -----------------------------
% Load Necessary Parameters
% -----------------------------

height = 100;  % Change this to the desired stack height
num_years = 5; % Set the number of years to simulate

% Define wind data filenames for each year
wind_files = cell(num_years, 1);
for year = 1:num_years
    wind_files{year} = sprintf('wind_results_%dm_combined_%dyear.mat', height, year);
end

% Initialize variables to store Nt and delta_t for consistency checks
Nt_all = zeros(num_years,1);
delta_t_all = zeros(num_years,1);

% Load grid data from the Excel file (assuming it's the same for all years)
data = readtable('RSudeste_Ventos.xlsx');

% Extract grid information from data
lat = data.LAT;             % Latitudes in degrees
lon = data.LONG;            % Longitudes in degrees
GridX = 189;                % Number of grid points in X-direction
GridY = 127;                % Number of grid points in Y-direction

% Earth's mean radius
R = 6371000;              % Radius in meters

% Geographic limits
lat_min = min(lat);         % Minimum latitude in degrees
lat_max = max(lat);         % Maximum latitude in degrees
lon_min = min(lon);         % Minimum longitude in degrees
lon_max = max(lon);         % Maximum longitude in degrees

% Calculate the average latitude for proper scaling
lat_avg = (lat_max + lat_min) / 2;  % Average latitude in degrees

% Convert longitude difference to meters (delta_x)
delta_x = ((lon_max - lon_min) * (pi / 180) * R * cosd(lat_avg)) / (GridX - 1);  % Grid spacing in X-direction (meters)

% Convert latitude difference to meters (delta_y)
delta_y = ((lat_max - lat_min) * (pi / 180) * R) / (GridY - 1);  % Grid spacing in Y-direction (meters)

% -----------------------------
% CFD Parameters
% -----------------------------

Dm = 0;                     % Molecular diffusivity is negligible
k_dep = 0.001 / height;     % Deposition rate constant (s⁻¹)
rho_air = 1.225;            % Air density (kg/m³)

% Calculate the emission rate
M_total = 3277267; % Total annual emission in kg/year
M_rate = M_total / (365 * 24 * 3600); % Mass emission rate in kg/s
m_air = rho_air * delta_x * delta_y * height; % Mass of air in grid cell in kg
source_emission_rate = M_rate / m_air; % Source emission rate (kg pollutant/kg air/s)

% -----------------------------
% Initialize Matrices
% -----------------------------

C_current = zeros(GridX, GridY);        % Pollutant concentration at current time step
deposition_total = zeros(GridX, GridY); % Total deposition over all years
concentration_total = zeros(GridX, GridY);  % Total concentration over all years

% For calculating dC/dM
total_time_steps = 0; % Will accumulate based on Nt of each year

% Define Visualization Times
num_steps_per_hour = floor(3600 / 1);  % Assuming delta_t = 1 second initially
visualization_times = [];  % Will be updated per year

% Precompute Mapping from Data Indices to Grid Indices
% Assuming data.X and data.Y are provided and correspond to grid points
% Adjust as necessary based on actual data structure
lin_idx = sub2ind([GridX, GridY], data.X, data.Y);

% -----------------------------
% Initialize Combined Results
% -----------------------------

% Initialize combined totals
deposition_total = zeros(GridX, GridY);
concentration_total_combined = zeros(GridX, GridY);

% -----------------------------
% CFD Simulation Loop
% -----------------------------

fprintf('Starting CFD simulation for %d years...\n', num_years);

for year = 1:num_years
    fprintf('Loading wind data for Year %d...\n', year);
    % Load the wind data for the current year using matfile
    mf = matfile(wind_files{year});
    
    % Validate consistency of Nt and delta_t
    Nt_all(year) = mf.Nt;                 % Total number of time steps per year
    delta_t_all(year) = mf.delta_tps;     % Time step size in seconds
    
    if year > 1
        if Nt_all(year) ~= Nt_all(1) || delta_t_all(year) ~= delta_t_all(1)
            error('Inconsistent Nt or delta_t across wind files.');
        end
    end
    
    Nt = Nt_all(year);
    delta_t = delta_t_all(year);
    
    % Update total_time_steps
    total_time_steps = total_time_steps + Nt;
    
    % Update num_steps_per_hour based on delta_t
    num_steps_per_hour = floor(3600 / delta_t);  % Number of time steps per hour
    
    % Calculate number of hours per year
    num_hours_per_year = floor(Nt / num_steps_per_hour);
    
    % Update visualization_times
    if year == 1
        visualization_times = num_steps_per_hour : num_steps_per_hour : Nt * num_years;
    else
        visualization_times = [visualization_times, total_time_steps - Nt + num_steps_per_hour : num_steps_per_hour : total_time_steps];
    end
    
    % Initialize yearly totals
    deposition_total_yearly = zeros(GridX, GridY);
    concentration_total_yearly = zeros(GridX, GridY);
    
    % Initialize hourly snapshots for the new year
    deposition_hourly = zeros(GridX, GridY);  % Hourly deposition accumulation
    deposition_hourly_snapshots = cell(num_hours_per_year, 1);
    concentration_hourly_snapshots = cell(num_hours_per_year, 1);
    hour_idx = 0;
    
    fprintf('Simulating Year %d...\n', year);
    
    for t_in_year = 1:Nt
        t = (year-1)*Nt + t_in_year; % Global time step
        
        % Display progress every 10000 time steps
        if mod(t, 10000) == 0
            fprintf('Processing Year %d: Time step %d of %d\n', year, t_in_year, Nt);
        end
        
        % Read wind data for current time step
        u_wind_t = mf.u_wind(:, t_in_year);    % Wind velocity in X-direction
        v_wind_t = mf.v_wind(:, t_in_year);    % Wind velocity in Y-direction
        
        % Map wind data to grid
        u_wind_grid = zeros(GridX, GridY);
        v_wind_grid = zeros(GridX, GridY);
        u_wind_grid(lin_idx) = u_wind_t;
        v_wind_grid(lin_idx) = v_wind_t;
        
        % Calculate wind speed magnitude and turbulent diffusivity Dt
        V_mag_grid = sqrt(u_wind_grid.^2 + v_wind_grid.^2);
        Dt_grid = V_mag_grid * height;
        
        % Initialize C_next
        C_next = zeros(GridX, GridY);
        
        % CFD calculations
        for i = 2 : GridX - 1
            for j = 2 : GridY - 1
                % Wind speed magnitude and turbulent diffusivity at current point
                V_mag = V_mag_grid(i, j);
                Dt = Dt_grid(i, j);
                
                % Advection Terms
                % X-direction (upwind differencing)
                if u_wind_grid(i, j) >= 0
                    derx = u_wind_grid(i, j) * (C_current(i, j) - C_current(i - 1, j)) / delta_x;
                else
                    derx = u_wind_grid(i, j) * (C_current(i + 1, j) - C_current(i, j)) / delta_x;
                end
                
                % Y-direction (upwind differencing)
                if v_wind_grid(i, j) >= 0
                    dery = v_wind_grid(i, j) * (C_current(i, j) - C_current(i, j - 1)) / delta_y;
                else
                    dery = v_wind_grid(i, j) * (C_current(i, j + 1) - C_current(i, j)) / delta_y;
                end
                
                % Diffusion Term (Laplacian)
                laplacian = (C_current(i + 1, j) + C_current(i - 1, j) - 2 * C_current(i, j)) / delta_x^2 + ...
                            (C_current(i, j + 1) + C_current(i, j - 1) - 2 * C_current(i, j)) / delta_y^2;
                
                % Source Term
                if i == 95 && j == 65
                    source_term = source_emission_rate;
                else
                    source_term = 0;
                end
                
                % Update Concentration
                C_next(i, j) = C_current(i, j) + delta_t * ( ...
                               - (derx + dery) + Dt * laplacian + source_term );
                
                % Deposition Calculation based on C_current
                delta_C_dep = C_current(i, j) * k_dep * delta_t;
                C_next(i, j) = C_next(i, j) - delta_C_dep;
                deposition_sec = delta_C_dep * rho_air;
                
                % Ensure non-negative concentrations
                if C_next(i, j) < 0
                    C_next(i, j) = 0;
                end
                
                % Update depositions
                deposition_total_yearly(i, j) = deposition_total_yearly(i, j) + deposition_sec;
                deposition_total(i, j) = deposition_total(i, j) + deposition_sec;
                deposition_hourly(i, j) = deposition_hourly(i, j) + deposition_sec;
            end
        end
        
        % Boundary Conditions (Neumann)
        C_next(1, :) = C_next(2, :);
        C_next(GridX, :) = C_next(GridX - 1, :);
        C_next(:, 1) = C_next(:, 2);
        C_next(:, GridY) = C_next(:, GridY - 1);
        
        % Prepare for Next Time Step
        C_current = C_next;
        
        % Accumulate total concentration
        concentration_total_yearly = concentration_total_yearly + C_current;
        concentration_total_combined = concentration_total_combined + C_current;
        
        % Save Snapshots at Visualization Times (Hourly)
        if mod(t_in_year, num_steps_per_hour) == 0
            hour_idx = hour_idx + 1;
            deposition_hourly_snapshots{hour_idx} = deposition_hourly;
            deposition_hourly = zeros(GridX, GridY);  % Reset hourly deposition

            concentration_hourly_snapshots{hour_idx} = C_current;  % Store current concentration snapshot
        end
    end
    
    % At the end of the year, save data
    fprintf('Year %d completed.\n', year);
    save_filename = sprintf('2VCCcfd_results_year%d_%dm.mat', year, height);
    save(save_filename, 'deposition_total_yearly', 'deposition_hourly_snapshots', ...
         'concentration_total_yearly', 'concentration_hourly_snapshots', ...
         'delta_x', 'delta_y', 'visualization_times', 'delta_t', 'lat_min', 'lat_max', ...
         'lon_min', 'lon_max', 'GridX', 'GridY', '-v7.3');
end

fprintf('CFD simulation for all %d years completed.\n', num_years);

% Compute the average concentration over all time steps and years
average_concentration = concentration_total_combined / total_time_steps;
average_concentration_kgm3 = average_concentration * rho_air; % kg/m³

% Calculate dC/dM as the average concentration per unit emission rate
dC_dM = average_concentration / M_rate; % Units: m³/kg

% Save the cumulative results over all years, including dC/dM
save_filename = sprintf('2VCCcfd_results_%dm_combined_%dyears.mat', height, num_years);
save(save_filename, 'Nt', 'deposition_total', 'average_concentration_kgm3', 'dC_dM', ...
     'delta_x', 'delta_y', 'delta_t', 'lat_min', 'lat_max', ...
     'lon_min', 'lon_max', 'GridX', 'GridY', '-v7.3');
