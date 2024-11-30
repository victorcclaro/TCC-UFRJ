clear all; close all; clc;

% Load Necessary Parameters
height = 100;  % Change this to the desired height (30, 50, 80, 100, 120, 150, 200)
year = 1;

% Visualization Loop
cfd_data_filename = sprintf('cfd_results_year%d_%dm.mat', year, height);
	
% Load the CFD results including snapshots
load(cfd_data_filename, 'concentration_hourly_snapshots', ...
    'visualization_times', 'delta_t', 'lat_min', 'lat_max', ...
    'lon_min', 'lon_max', 'GridX', 'GridY');

%% Retrieve the Map Image from WMS Service
% Find the OpenStreetMap layer
osmLayer = wmsfind('osm', 'SearchField', 'serverurl');
osmLayer = refine(osmLayer, 'osm');
layer = osmLayer(1);

% Define the image size
imageWidth = 2000;  
imageHeight = 1600;

% Retrieve the map image with specified latitude and longitude limits
[A, R] = wmsread(layer, 'Latlim', [lat_min lat_max], 'Lonlim', [lon_min lon_max], ...
    'ImageHeight', imageHeight, 'ImageWidth', imageWidth);

% Add Distance Scale/Legend
% Calculate the center latitude
center_lat = (lat_min + lat_max) / 2;

% Calculate degrees per kilometer
degrees_per_km_lat = 1 / 111.32;  % Degrees latitude per km
degrees_per_km_lon = 1 / (111.32 * cosd(center_lat));  % Degrees longitude per km

% Define the scale bar length in kilometers
scale_length_km = 100;  % 100 km scale bar

% Calculate the scale bar length in degrees of longitude
scale_bar_length_deg = scale_length_km * degrees_per_km_lon;

% Define offsets from the map edges (as a percentage of map dimensions)
offset_x_pct = 0.05;  % 5% offset from the right edge
offset_y_pct = 0.05;  % 5% offset from the bottom edge

% Calculate the offsets in degrees
offset_x_deg = (lon_max - lon_min) * offset_x_pct;
offset_y_deg = (lat_max - lat_min) * offset_y_pct;

% Define the starting point of the scale bar
x_start = lon_max - offset_x_deg - scale_bar_length_deg;
y_start = lat_min + offset_y_deg;

% Define the end point of the scale bar
x_end = x_start + scale_bar_length_deg;
y_end = y_start;

% Create the Figure and Display the Background Map
figure;

% Display the map with geographic referencing
hMap = imshow(A, 'XData', [lon_min lon_max], 'YData', [lat_min lat_max]);  
set(hMap, 'AlphaData', 0.5);  % Adjust transparency of the map

% Ensure that the aspect ratio is correct
axis on;        % Turn on axis
axis equal;     % Keep the geographic aspect ratio

hold on;        % Hold on to overlay the contour plot

% Plot the scale bar
plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', 2);
plot([x_start, x_start], [y_start - offset_y_deg * 0.1, y_start + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
plot([x_end, x_end], [y_end - offset_y_deg * 0.1, y_end + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
text(x_start + scale_bar_length_deg/2, y_start - offset_y_deg, ...
    [num2str(scale_length_km) ' km'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');

% Generate Meshgrid for Contour Plot
[X, Y] = meshgrid(linspace(lon_min, lon_max, GridX), linspace(lat_min, lat_max, GridY));

% Create the Contour Plot Using the First Snapshot
% Transpose the data to match the meshgrid orientation
[c, cont] = contourf(X, Y, concentration_hourly_snapshots{1}' * 1e9, 'LineColor', 'none');

% Set axis limits to match the geographic limits
set(gca, 'XLim', [lon_min lon_max], 'YLim', [lat_min lat_max]);

% Add a colorbar to show the value scale for the contour plot
colorbar;

% Set titles and labels
title('PM10 Concentration at Hour 1 [µg/m³]');
xlabel('Longitude');
ylabel('Latitude');

% Time Loop to Update the Contour Plot with New Data
for idx = 1:length(visualization_times)
    t = visualization_times(idx);
    
    % Calculate time in hours
    time_in_hours = t * delta_t / 3600;
    
    % Update the title with the current hour
    title(sprintf('PM10 Concentration at Hour %.0f [µg/m³]', time_in_hours));
    
    % Update the ZData of the contour plot for the current time step
    set(cont, 'ZData', concentration_hourly_snapshots{idx}' * 1e9);
    
    % Bring the contour plot to the bottom
    uistack(cont, 'bottom');
    
    % Pause to create an animation effect
    pause(0.001);  % Adjust the pause duration as needed
end

