clear all; close all; clc;

% Load data from Excel file
data = readtable('RSudeste_Ventos.xlsx');

% Extract relevant columns from the table
lat = data.LAT;
lon = data.LONG;

% Define compass points and their corresponding angles in degrees
compass_points = {'N', 'NNE', 'NE', 'NEE', 'E', 'ESE', ...
                  'SE', 'SSE', 'S', 'SSO', 'SO', 'OSO', ...
                  'O', 'ONO', 'NO', 'NNO'};

angles_deg = [0, 22.5, 45, 67.5, 90, 112.5, ...
              135, 157.5, 180, 202.5, 225, 247.5, ...
              270, 292.5, 315, 337.5];

% Wind direction frequencies for each compass point
wind_dir_freqs = [data.N, data.NNE, data.NE, data.NEE, data.E, data.ESE, ...
                  data.SE, data.SSE, data.S, data.SSO, data.SO, data.OSO, ...
                  data.O, data.ONO, data.NO, data.NNO];  % [num_locations x 16]

% Wind speed factors for each hour (F1h to F24h)
wind_factors = [data.F1h, data.F2h, data.F3h, data.F4h, data.F5h, data.F6h, ...
                data.F7h, data.F8h, data.F9h, data.F10h, data.F11h, data.F12h, ...
                data.F13h, data.F14h, data.F15h, data.F16h, data.F17h, data.F18h, ...
                data.F19h, data.F20h, data.F21h, data.F22h, data.F23h, data.F24h];

% Wind speed data at 100m only
wind_speeds = data.V_100m;

% Preallocate matrices to store results
num_locations = size(data, 1);  % Number of locations
num_hours = 24;  % Number of hours in the day
num_days = 365;  % Number of days

relative_wind = nan(num_locations, num_hours, 1);  % [num_locations x 24 x 1]

% Normalize wind direction frequencies for each location
freq_sums = sum(wind_dir_freqs, 2);  % Sum across compass points for each location
normalized_freqs = wind_dir_freqs ./ freq_sums;  % [num_locations x 16]

% Replace NaNs and zeros in normalized frequencies with equal probabilities
normalized_freqs(isnan(normalized_freqs)) = 1/16;
normalized_freqs(freq_sums == 0, :) = 1/16;

% Precompute angular distances between compass points
num_states = length(angles_deg);  % Should be 16
angular_distances = zeros(num_states, num_states);

for s = 1:num_states
    for s_next = 1:num_states
        diff = abs(angles_deg(s) - angles_deg(s_next));
        angular_distances(s, s_next) = min(diff, 360 - diff);
    end
end

% Define beta parameter (controls the strength of temporal correlation)
beta = 0.1;  % Adjust beta as needed to control temporal correlation

% Initialize cell array to store transition matrices
transition_matrices = cell(num_locations, 1);

% Construct transition matrices for each location
for loc = 1:num_locations
    T = zeros(num_states, num_states);
    for s = 1:num_states
        % Compute unnormalized transition probabilities
        T(s, :) = normalized_freqs(loc, :) .* exp(-beta * angular_distances(s, :));
        % Normalize to sum to 1
        T(s, :) = T(s, :) / sum(T(s, :));
    end
    % Store transition matrix for location
    transition_matrices{loc} = T;
end

% Loop over each location to process the data
for i = 1:num_locations
    % Average the wind speed factors for the 24 hours
    avg_wind_factor = mean(wind_factors(i, :), 'omitnan');

    % Calculate relative wind differences for each hour
    relative_factors = wind_factors(i, :) / avg_wind_factor;

    % Multiply the relative factors by the wind speed at 100m
    relative_wind(i, :, 1) = relative_factors .* wind_speeds(i);
end

% Time parameters
delta_tps = 300;  % Time step size in seconds (5 minutes)
num_steps_per_hour = 3600 / delta_tps;  % Number of time steps per hour (12)
Nt = num_days * num_hours * num_steps_per_hour;  % Total number of time steps

% Prepare to save wind data for 100m only
heights = 100;  % Heights in meters

% Define batch size (adjust based on available memory)
batch_size = 1000;  % Adjust this number based on your system's available memory
num_batches = ceil(num_locations / batch_size);

for h_idx = 1:length(heights)  % This will run only once
    height = heights;  % 100m
    fprintf('Processing wind data for height %dm...\n', height);

    % Initialize matfile object for combined data
    save_filename = sprintf('5wind_data_%dm_combined.mat', height);
    mf = matfile(save_filename, 'Writable', true);

    % Optionally, initialize arrays in the MAT-file
    mf.u_wind = [];  % Initialize as empty
    mf.v_wind = [];
    mf.Nt = Nt;
    mf.delta_tps = delta_tps;

    % Generate wind directions for each hour over all days
    num_total_hours = num_hours * num_days;

    for batch = 1:num_batches
        idx_start = (batch - 1) * batch_size + 1;
        idx_end = min(batch * batch_size, num_locations);
        batch_indices = idx_start:idx_end;
        num_locations_in_batch = length(batch_indices);

        fprintf('  Processing batch %d of %d (%d locations)...\n', batch, num_batches, num_locations_in_batch);

        % Initialize arrays for this batch
        wind_directions_per_hour = zeros(num_locations_in_batch, num_total_hours);  % [batch_size x num_total_hours]

        % Generate wind directions for each location in the batch
        for loc_idx = 1:num_locations_in_batch
            loc = batch_indices(loc_idx);
            T = transition_matrices{loc};
            dir_idx = randsample(1:num_states, 1, true, normalized_freqs(loc, :));
            wind_directions_per_hour(loc_idx, 1) = angles_deg(dir_idx);

            for h_idx_time = 2:num_total_hours
                current_state = dir_idx;
                transition_probs = T(current_state, :);
                dir_idx = randsample(1:num_states, 1, true, transition_probs);
                wind_directions_per_hour(loc_idx, h_idx_time) = angles_deg(dir_idx);
            end
        end

        % Initialize simulated wind arrays for the batch
        simulated_wind = nan(num_locations_in_batch, Nt);
        simulated_wind_dir = nan(num_locations_in_batch, Nt);

        % Populate simulated wind data for the batch
        for day = 1:num_days
            for h = 1:num_hours
                t_start = (day - 1) * num_hours * num_steps_per_hour + (h - 1) * num_steps_per_hour + 1;
                t_end = t_start + num_steps_per_hour - 1;

                wind_speed_at_height = relative_wind(batch_indices, h, h_idx);  % [num_locations_in_batch x 1]
                hour_idx = (day - 1) * num_hours + h;
                wind_dir_at_hour = wind_directions_per_hour(:, hour_idx);  % [num_locations_in_batch x 1]

                simulated_wind(:, t_start:t_end) = repmat(wind_speed_at_height, 1, num_steps_per_hour);
                simulated_wind_dir(:, t_start:t_end) = repmat(wind_dir_at_hour, 1, num_steps_per_hour);
            end
        end

        % Calculate wind vector components (u, v) based on sampled wind directions
        u_wind = simulated_wind .* cosd(simulated_wind_dir);  % [num_locations_in_batch x Nt]
        v_wind = simulated_wind .* sind(simulated_wind_dir);  % [num_locations_in_batch x Nt]

        % Write batch data into the matfile at the correct indices
        mf.u_wind(batch_indices, 1:Nt) = u_wind;
        mf.v_wind(batch_indices, 1:Nt) = v_wind;

        % Clear variables to free up memory
        clear u_wind v_wind simulated_wind simulated_wind_dir wind_directions_per_hour;
    end

    % Save additional variables if needed
    mf.Nt = Nt;
    mf.delta_tps = delta_tps;

    fprintf('Wind data processing and saving completed for height %dm.\n', height);
end

fprintf('All wind data processing and saving completed.\n');
