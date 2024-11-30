clear all; close all; clc;

% Define the grid and boundaries
% Read grid data from the Excel file
data = readtable('RSudeste_Ventos.xlsx');

% Extract grid information from data
lat = data.LAT;             % Latitudes in degrees
lon = data.LONG;            % Longitudes in degrees

% Geographic limits
lat_min = min(lat);         % Minimum latitude in degrees
lat_max = max(lat);         % Maximum latitude in degrees
lon_min = min(lon);         % Minimum longitude in degrees
lon_max = max(lon);

% Create a grid with 189 x 127 points
num_lat_points = 127;
num_lon_points = 189;


% Generate latitude and longitude vectors for the grid
lat_grid = linspace(lat_min, lat_max, num_lat_points);
lon_grid = linspace(lon_min, lon_max, num_lon_points);

% Generate a meshgrid for the latitudes and longitudes
[lon_mesh, lat_mesh] = meshgrid(lon_grid, lat_grid);

%% Full Map with Highlighted Cell (95, 64)

% Create a figure for the full map
figure('Position', [100, 100, 1200, 1000]);

% Set up map projection using Mapping Toolbox
axesm('mercator', 'MapLatLimit', [lat_min lat_max], 'MapLonLimit', [lon_min lon_max]);
axis off;

% Find an OpenStreetMap WMS layer to use for the background map
osmLayer = wmsfind('osm', 'SearchField', 'serverurl');  
osmLayer = refine(osmLayer, 'osm');
layer = osmLayer(1); 

% Retrieve the map data from the WMS server
imageWidth = 2000;  
imageHeight = 1600;  

[A, R] = wmsread(layer, 'Latlim', [lat_min lat_max], 'Lonlim', [lon_min lon_max], ...
    'ImageHeight', imageHeight, 'ImageWidth', imageWidth);

% Display the high-resolution OpenStreetMap as the background
geoshow(A, R);
hold on;

% Overlay the black grid on the map
plotm(lat_mesh, lon_mesh, 'k-', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]);  
plotm(lat_mesh', lon_mesh', 'k-', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]);  

% Highlight only the cell at position (95, 64) in red
cell_lon_start = lon_grid(94);
cell_lon_end = lon_grid(95);
cell_lat_start = lat_grid(63);
cell_lat_end = lat_grid(64);

% Plot the red boundary for the specific cell (95, 64)
plotm([cell_lat_start, cell_lat_start], [cell_lon_start, cell_lon_end], 'r-', 'LineWidth', 1.5);
plotm([cell_lat_end, cell_lat_end], [cell_lon_start, cell_lon_end], 'r-', 'LineWidth', 1.5);
plotm([cell_lat_start, cell_lat_end], [cell_lon_start, cell_lon_start], 'r-', 'LineWidth', 1.5);
plotm([cell_lat_start, cell_lat_end], [cell_lon_end, cell_lon_end], 'r-', 'LineWidth', 1.5);

% Calculate and plot the scale bar for the full map
center_lat = (lat_min + lat_max) / 2;
degrees_per_km_lon = 1 / (111.32 * cosd(center_lat));
scale_length_km = 100;
scale_bar_length_deg = scale_length_km * degrees_per_km_lon;
offset_x_pct = 0.05;  
offset_y_pct = 0.05;  

x_start = lon_max - (lon_max - lon_min) * offset_x_pct - scale_bar_length_deg;
y_start = lat_min + (lat_max - lat_min) * offset_y_pct;
x_end = x_start + scale_bar_length_deg;
y_end = y_start;

% Plot the scale bar
plotm([y_start, y_end], [x_start, x_end], 'k-', 'LineWidth', 2);
plotm([y_start - (lat_max - lat_min) * 0.005, y_start + (lat_max - lat_min) * 0.005], [x_start, x_start], 'k-', 'LineWidth', 2);
plotm([y_end - (lat_max - lat_min) * 0.005, y_end + (lat_max - lat_min) * 0.005], [x_end, x_end], 'k-', 'LineWidth', 2);
textm(y_start - (lat_max - lat_min) * 0.02, (x_start + x_end) / 2, [num2str(scale_length_km) ' km'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');

% Invert y-axis labels (Latitude) every 10 grid points on the left edge
for j = 1:10:num_lat_points
    % Left edge labels with specific grid numbers
    inverted_lat_index = num_lat_points - j + 1; % Adjust for inversion
    textm(lat_grid(inverted_lat_index), lon_min, sprintf('%d', inverted_lat_index), ...
          'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
          'FontSize', 8, 'FontWeight', 'bold');
end

% Add X-axis labels (Longitude) every 10 grid points only on the bottom edge
for i = 1:10:num_lon_points
    % Bottom edge labels with specific grid numbers
    textm(lat_min, lon_grid(i), sprintf('%d', i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'FontWeight', 'bold');
end
% Save the zoomed map with scale bar
print('map_with_highlighted_cells_and_scale_bar', '-dpng', '-r600');  


%% Zoomed Map on Region (90:100, 59:69)

% Define the zoomed region limits
zoom_lat_min = lat_grid(59);
zoom_lat_max = lat_grid(69);
zoom_lon_min = lon_grid(90);
zoom_lon_max = lon_grid(100);

% Create a new figure for the zoomed map
figure('Position', [100, 100, 1200, 1200]);

% Set up map projection for the zoomed region
axesm('mercator', 'MapLatLimit', [zoom_lat_min zoom_lat_max], 'MapLonLimit', [zoom_lon_min zoom_lon_max]);
axis off;

% Retrieve zoomed map data from the WMS server
[A_zoom, R_zoom] = wmsread(layer, 'Latlim', [zoom_lat_min zoom_lat_max], 'Lonlim', [zoom_lon_min zoom_lon_max], ...
    'ImageHeight', imageHeight, 'ImageWidth', imageWidth);

% Display the zoomed-in map
geoshow(A_zoom, R_zoom);
hold on;

% Plot the grid on the zoomed map
[lon_mesh_zoom, lat_mesh_zoom] = meshgrid(linspace(zoom_lon_min, zoom_lon_max, 11), linspace(zoom_lat_min, zoom_lat_max, 11));
plotm(lat_mesh_zoom, lon_mesh_zoom, 'k-', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]);
plotm(lat_mesh_zoom', lon_mesh_zoom', 'k-', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]);

% Define the latitude and longitude boundaries for the highlighted cells
cell_coords = {
    [96, 64], 'b';   % Cell (96, 64) in blue
    [95, 65], 'b';   % Cell (95, 65) in blue
    [96, 65], 'b'    % Cell (96, 65) in blue
    [95, 64], 'r';   % Cell (95, 64) in red
};

% Highlight specific cells in red or blue within the zoomed region
for k = 1:size(cell_coords, 1)
    x_idx = cell_coords{k, 1}(1);
    y_idx = cell_coords{k, 1}(2);
    color = cell_coords{k, 2};

    if x_idx >= 90 && x_idx <= 100 && y_idx >= 59 && y_idx <= 69
        cell_lon_start = lon_grid(x_idx - 1);
        cell_lon_end = lon_grid(x_idx);
        cell_lat_start = lat_grid(y_idx - 1);
        cell_lat_end = lat_grid(y_idx);

        % Plot cell boundaries in specified color
        plotm([cell_lat_start, cell_lat_start], [cell_lon_start, cell_lon_end], color, 'LineWidth', 1.5);
        plotm([cell_lat_end, cell_lat_end], [cell_lon_start, cell_lon_end], color, 'LineWidth', 1.5);
        plotm([cell_lat_start, cell_lat_end], [cell_lon_start, cell_lon_start], color, 'LineWidth', 1.5);
        plotm([cell_lat_start, cell_lat_end], [cell_lon_end, cell_lon_end], color, 'LineWidth', 1.5);
    end
end

% Scale bar for the zoomed map
zoom_center_lat = (zoom_lat_min + zoom_lat_max) / 2;
zoom_degrees_per_km_lon = 1 / (111.32 * cosd(zoom_center_lat));
zoom_scale_length_km = 10;  
zoom_scale_bar_length_deg = zoom_scale_length_km * zoom_degrees_per_km_lon;

x_start_zoom = zoom_lon_max - (zoom_lon_max - zoom_lon_min) * offset_x_pct - zoom_scale_bar_length_deg;
y_start_zoom = zoom_lat_min + (zoom_lat_max - zoom_lat_min) * offset_y_pct;
x_end_zoom = x_start_zoom + zoom_scale_bar_length_deg;
y_end_zoom = y_start_zoom;

% Plot the scale bar
plotm([y_start_zoom, y_end_zoom], [x_start_zoom, x_end_zoom], 'k-', 'LineWidth', 2);
plotm([y_start_zoom - (zoom_lat_max - zoom_lat_min) * 0.005, y_start_zoom + (zoom_lat_max - zoom_lat_min) * 0.005], [x_start_zoom, x_start_zoom], 'k-', 'LineWidth', 2);
plotm([y_end_zoom - (zoom_lat_max - zoom_lat_min) * 0.005, y_end_zoom + (zoom_lat_max - zoom_lat_min) * 0.005], [x_end_zoom, x_end_zoom], 'k-', 'LineWidth', 2);
textm(y_start_zoom - (zoom_lat_max - zoom_lat_min) * 0.02, (x_start_zoom + x_end_zoom) / 2, [num2str(zoom_scale_length_km) ' km'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');

% Calculate the spacing between points in the zoomed region
lat_spacing = (zoom_lat_max - zoom_lat_min) / (length(lat_mesh_zoom) - 1);
lon_spacing = (zoom_lon_max - zoom_lon_min) / (length(lon_mesh_zoom) - 1);

% Add inverted Y-axis labels (Latitude) every 10 grid points on the left edge, moving them down
for j = 1:1:num_lat_points
    inverted_lat_index = num_lat_points - j + 1; % Adjust for inversion

    if lat_grid(inverted_lat_index) >= zoom_lat_min && lat_grid(inverted_lat_index) <= zoom_lat_max
        % Move labels down by subtracting spacing
        adjusted_lat = lat_grid(inverted_lat_index) - lat_spacing / 2;
        
        % Use inverted_lat_index for the label value to reflect the inverted order
        textm(adjusted_lat, zoom_lon_min, sprintf('%d', inverted_lat_index), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 8, 'FontWeight', 'bold');
    end
end


% Add X-axis labels (Longitude) every 10 grid points only on the bottom edge, moving them left
for i = 1:1:num_lon_points
    if lon_grid(i) >= zoom_lon_min && lon_grid(i) <= zoom_lon_max
        % Move labels left by subtracting spacing
        adjusted_lon = lon_grid(i) - lon_spacing / 2;
        textm(zoom_lat_min, adjusted_lon, sprintf('%d', i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8, 'FontWeight', 'bold');
    end
end

hold off;

% Save the zoomed map with scale bar
print('zoomed_map_with_highlighted_cells_and_scale_bar', '-dpng', '-r600');  
