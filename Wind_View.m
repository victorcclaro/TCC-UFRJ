clear all; close all; clc;

% Load necessary parameters from the wind simulation data
height = 100;  % Change this to the desired height (30, 50, 80, 100, 120, 150, 200)

% Correct wind data filename
wind_data_filename = sprintf('wind_results_%dm_combined_1year.mat', height);
fprintf('Loading wind data from %s...\n', wind_data_filename);

% Use matfile to access data incrementally
mf = matfile(wind_data_filename);

% Load small variables
Nt = mf.Nt;
delta_tps = mf.delta_tps;

% Read 'data' directly from the Excel file
data = readtable('RSudeste_Ventos.xlsx');

% Extract grid information from data
lat = data.LAT;
lon = data.LONG;
GridX = 189;  % Number of cells in X-direction
GridY = 127;  % Number of cells in Y-direction

% Precompute mapping from data indices to grid indices
% Ensure 'data.X' and 'data.Y' exist or compute them
if ~ismember('X', data.Properties.VariableNames) || ~ismember('Y', data.Properties.VariableNames)
    % Compute 'data.X' and 'data.Y' based on lat/lon
    lat_min = min(lat);
    lat_max = max(lat);
    lon_min = min(lon);
    lon_max = max(lon);

    data.X = round((data.LONG - lon_min) / (lon_max - lon_min) * (GridX - 1)) + 1;
    data.Y = round((data.LAT - lat_min) / (lat_max - lat_min) * (GridY - 1)) + 1;

    % Ensure indices are within grid bounds
    data.X = min(max(data.X, 1), GridX);
    data.Y = min(max(data.Y, 1), GridY);
end

% Precompute linear indices for mapping
lin_idx = sub2ind([GridX, GridY], data.X, data.Y);

% Define visualization times
num_steps_per_hour = 3600 / delta_tps;  % Number of time steps per hour
visualization_times = num_steps_per_hour : num_steps_per_hour : Nt;

% Retrieve the map image from WMS service with specified lat/lon limits
lat_min = min(lat);
lat_max = max(lat);
lon_min = min(lon);
lon_max = max(lon);

% Retrieve the map image from WMS service with specified lat/lon limits
osmLayer = wmsfind('osm', 'SearchField', 'serverurl');
osmLayer = refine(osmLayer, 'osm');
layer = osmLayer(1);

imageWidth = 2000;  
imageHeight = 1600;

% Retrieve the map image from WMS service
[A, R] = wmsread(layer, 'Latlim', [lat_min lat_max], 'Lonlim', [lon_min lon_max], ...
    'ImageHeight', imageHeight, 'ImageWidth', imageWidth);

% Generate a meshgrid for the contour plot
[X_grid, Y_grid] = meshgrid(linspace(lon_min, lon_max, GridX), linspace(lat_min, lat_max, GridY));

%% Visualization: Wind vector field animation over the defined hours
figure;

% Display the map with geographic referencing (using XData and YData)
hMap = imshow(A, 'XData', [lon_min lon_max], 'YData', [lat_min lat_max]);  
set(hMap, 'AlphaData', 0.5);  % Adjust transparency of the map
axis on;  
axis equal;  
hold on;

% Add the scale bar
% Calculate center latitude for accurate scaling
center_lat = (lat_min + lat_max) / 2;
degrees_per_km_lon = 1 / (111.32 * cosd(center_lat));

% Define the scale bar length in kilometers
scale_length_km = 100;  % 100 km scale bar

% Calculate the scale bar length in degrees of longitude
scale_bar_length_deg = scale_length_km * degrees_per_km_lon;

% Define offsets from the map edges for the scale bar
offset_x_pct = 0.05;  % 5% offset from the right edge
offset_y_pct = 0.05;  % 5% offset from the bottom edge

% Calculate the offsets in degrees
offset_x_deg = (lon_max - lon_min) * offset_x_pct;
offset_y_deg = (lat_max - lat_min) * offset_y_pct;

% Define the starting and ending points of the scale bar
x_start = lon_max - offset_x_deg - scale_bar_length_deg;
y_start = lat_max - offset_y_deg;
x_end = x_start + scale_bar_length_deg;
y_end = y_start;

% Plot the scale bar
plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', 2);
plot([x_start, x_start], [y_start - offset_y_deg * 0.1, y_start + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
plot([x_end, x_end], [y_end - offset_y_deg * 0.1, y_end + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
text((x_start + x_end)/2, y_start - offset_y_deg * 0.5, ...
    [num2str(scale_length_km) ' km'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');

% Plot the city location
Xcity = 95;  % Grid index
Ycity = 64;  % Grid index
lon_city = lon_min + (Xcity - 1) / (GridX - 1) * (lon_max - lon_min);
lat_city = lat_min + (Ycity - 1) / (GridY - 1) * (lat_max - lat_min);
plot(lon_city, lat_city, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');

% Set axis limits to match the geographic limits
set(gca, 'XLim', [lon_min lon_max], 'YLim', [lat_min lat_max]);

for hour = 1:length(visualization_times)
    idx = visualization_times(hour);  % Index for the current time in the simulation

    % Read wind data for the current time step
    u_wind_t = mf.u_wind(:, idx);  % size [num_locations x 1]
    v_wind_t = mf.v_wind(:, idx);

    % Map wind data to grid
    uw = zeros(GridX, GridY);
    vw = zeros(GridX, GridY);

    % Assign wind data to grid locations
    uw(lin_idx) = u_wind_t;
    vw(lin_idx) = v_wind_t;

    % Plot wind vectors (quiver) on top of the map
    hQuiver = quiver(X_grid, Y_grid, uw', vw', 'b', 'LineWidth', 0.1, 'AutoScaleFactor', 0.8, 'MaxHeadSize', 0.05);

    title(['Wind Vector Field (Height: ', num2str(height), 'm) - Hour ', num2str(hour)]);
    xlabel('Longitude');
    ylabel('Latitude');

    pause(0.1);  % Pause to create the animation effect

    % Clear the current quiver to prevent overlap in the next iteration
    delete(hQuiver);
end

hold off;
