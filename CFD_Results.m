clear all; close all; clc;

% -----------------------------
% 1. Define Parameters
% -----------------------------

height = 100;  % Change this to the desired height (30, 50, 80, 100, 120, 150, 200)
start_year = 1;
end_year = 5;

% Define emission site central indices
central_row = 95;  % Adjust as needed
central_col = 64;  % Adjust as needed

% Preallocate storage for all years
all_avg_conc = [];       % Concatenated average concentrations over all years
all_time_in_hours = [];  % Concatenated time vectors
yearly_averages = zeros(end_year, 1);  % Yearly average concentrations

% Initialize variables to store problematic timesteps across all years
all_problematic_timesteps = cell(end_year,1);

% -----------------------------
% 2. Process Each Year
% -----------------------------

for year = start_year:end_year
    fprintf('Processing Year %d...\n', year);
    
    % Construct the filename based on year and height
    cfd_data_filename = sprintf('cfd_results_year%d_%dm.mat', year, height);
    
    % Check if the file exists
    if ~isfile(cfd_data_filename)
        warning('File %s does not exist. Skipping Year %d.', cfd_data_filename, year);
        continue;
    end
    
    % Load the CFD results including snapshots and relevant parameters
    load(cfd_data_filename, 'concentration_hourly_snapshots', ...
        'visualization_times', 'delta_t', 'lat_min', 'lat_max', ...
        'lon_min', 'lon_max', 'GridX', 'GridY');
    
    % -----------------------------
    % 3. Define the 2x2 Grid
    % -----------------------------
    
    % Define the range for rows and columns to form a 2x2 grid
    rows = (central_row):(central_row+1);  % 95, 96
    cols = (central_col):(central_col+1);  % 64, 65
    
    % Ensure that the indices do not exceed the grid boundaries
    rows = max(1, rows); 
    rows = min(GridX, rows);
    
    cols = max(1, cols); 
    cols = min(GridY, cols);
    
    % -----------------------------
    % 4. Verify Snapshot Dimensions
    % -----------------------------
    
    % Initialize variables to store maximum rows and columns
    max_rows = 0;
    max_cols = 0;
    
    % Determine the number of timesteps for the current year
    num_timesteps = length(concentration_hourly_snapshots);
    
    % Loop through all snapshots to find maximum dimensions
    for idx = 1:num_timesteps
        current_snapshot = concentration_hourly_snapshots{idx};
        [rows_snapshot, cols_snapshot] = size(current_snapshot);
        
        if rows_snapshot > max_rows
            max_rows = rows_snapshot;
        end
        if cols_snapshot > max_cols
            max_cols = cols_snapshot;
        end
    end
    
    fprintf('Maximum snapshot dimensions for Year %d: %d rows x %d columns\n', year, max_rows, max_cols);
    fprintf('Defined grid indices for Year %d: rows %d-%d, cols %d-%d\n', year, min(rows), max(rows), min(cols), max(cols));
    
    % Optional: Compare with GridX and GridY
    if max_rows > GridX || max_cols > GridY
        warning('Maximum snapshot dimensions exceed GridX/GridY for Year %d.', year);
    end
    
    % -----------------------------
    % 5. Calculate the Averages with Enhanced Boundary Checks
    % -----------------------------
    
    % Preallocate an array to store average concentrations for the current year
    avg_conc = zeros(num_timesteps, 1);
    
    % Initialize an array to store problematic timesteps for the current year
    problematic_timesteps = [];
    
    % Loop through each timestep to compute the average concentration
    for idx = 1:num_timesteps
        % Extract the current snapshot and apply unit conversion
        current_snapshot = concentration_hourly_snapshots{idx} * 1e9;  % Adjust if necessary
        
        % Get the size of the current snapshot
        [snapshot_rows, snapshot_cols] = size(current_snapshot);
        
        % Adjust rows and cols for this snapshot
        valid_rows = rows(rows <= snapshot_rows);
        valid_rows = valid_rows(valid_rows >= 1);
        
        valid_cols = cols(cols <= snapshot_cols);
        valid_cols = valid_cols(valid_cols >= 1);
        
        % Check if the 2x2 grid is fully available
        if length(valid_rows) < 2 || length(valid_cols) < 2
            warning('Year %d, Timestep %d: 2x2 grid exceeds snapshot boundaries. Adjusting to available data.', year, idx);
            problematic_timesteps(end+1) = idx;
        end
        
        % Extract the 2x2 grid around the emission site
        grid_values = current_snapshot(valid_rows, valid_cols);
        
        % Check if grid_values is non-empty
        if isempty(grid_values)
            warning('Year %d, Timestep %d: No valid grid data available.', year, idx);
            avg_conc(idx) = NaN;  % Assign NaN to indicate missing data
        else
            % Calculate the average concentration over the grid
            avg_val = mean(grid_values, 'all');
            
            % Optionally, round the average value
            avg_conc(idx) = round(avg_val);
        end
    end
    
    % Store problematic timesteps for the current year
    all_problematic_timesteps{year} = problematic_timesteps;
    
    % -----------------------------
    % 6. Aggregate Data Across Years
    % -----------------------------
    
    % Calculate the time vector in hours for the current year
    time_in_hours = (0:num_timesteps-1)' + (year-1)*8760;  % Offset by previous years
    
    % Concatenate the data
    all_avg_conc = [all_avg_conc; avg_conc];
    all_time_in_hours = [all_time_in_hours; time_in_hours];
    
    % Calculate and store the yearly average for the current year
    yearly_avg = mean(avg_conc, 'omitnan');
    yearly_averages(year) = yearly_avg;
    
    fprintf('Year %d: Yearly Average Concentration = %.2f\n', year, yearly_avg);
end

% -----------------------------
% 7. Plot the Time Series Over All Years
%% -----------------------------

figure;
plot(all_time_in_hours, all_avg_conc, 'b-', 'LineWidth', 1.5);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Average Concentration (µg/m³)', 'FontSize', 14);
title('Average Concentration at Emission Area Over 10 Years', 'FontSize', 16);
grid on;
hold on;

% Add vertical lines and year labels at each year transition
for yr = 1:end_year
    transition_time = yr * 8760;  % Number of hours in a year
    xline(transition_time, '--k', 'LineWidth', 1.5);
    
    % Add year labels slightly above the maximum concentration
    text(transition_time, max(all_avg_conc, [], 'omitnan')*1.02, sprintf('Year %d', yr), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Rotation', 90, ...
        'VerticalAlignment', 'bottom');
end

% Adjust x-axis limits to encompass all years
xlim([0, end_year * 8760]);

% Adjust y-axis limits based on data
min_conc = min(all_avg_conc, [], 'omitnan');
max_conc = max(all_avg_conc, [], 'omitnan');
ylim([min_conc*0.95, max_conc*1.05]);  % Add some padding

hold off;

% -----------------------------
% 8. Plot the Yearly Averages as a Bar Chart
%% -----------------------------

figure;
bar(start_year:end_year, yearly_averages, 0.6, 'FaceColor', [0.2 0.6 0.5]);
xlabel('Year', 'FontSize', 14);
ylabel('Yearly Average Concentration (µg/m³)', 'FontSize', 14);
title('Yearly Average Concentration at Emission Area', 'FontSize', 16);
grid on;

% Add data labels on top of each bar
for yr = 1:end_year
    if yr <= length(yearly_averages)
        text(yr, yearly_averages(yr), sprintf('%.2f', yearly_averages(yr)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12, 'FontWeight', 'bold');
    end
end

% Adjust x-axis limits
xlim([0.5, end_year + 0.5]);

% -----------------------------
% 9. Optional: Save the Plots
% -----------------------------

% Save the time series plot
% saveas(gcf, 'average_concentration_over_10_years.png');

% Save the yearly averages bar chart
% saveas(gcf, 'yearly_average_concentration.png');
