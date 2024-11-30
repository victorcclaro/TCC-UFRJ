clear all; close all; clc;

%% 1. Read Grid Data from the Excel File
data = readtable('RSudeste_Ventos.xlsx');

% Extract grid information from data
lat = data.LAT;             % Latitudes in degrees
lon = data.LONG;            % Longitudes in degrees
GridX = 189;                % Number of grid points in X-direction
GridY = 127;                % Number of grid points in Y-direction

%% 2. Define Geographic Limits
lat_min = min(lat);         % Minimum latitude in degrees
lat_max = max(lat);         % Maximum latitude in degrees
lon_min = min(lon);         % Minimum longitude in degrees
lon_max = max(lon);         % Maximum longitude in degrees

%% 3. Read Population Data
population_data_filename = 'gpw_v4_population_count_rev11_2020_30_sec.tif';

% Read the population data and the referencing object
[pop_data_full, pop_R] = readgeoraster(population_data_filename, 'OutputType', 'double');

% Convert any missing or negative values to zero
pop_data_full(pop_data_full < 0) = 0;
pop_data_full(isnan(pop_data_full)) = 0;

% Display raster referencing information
disp('Population Data Spatial Referencing Information:');
disp(pop_R);

%% 4. Define Latitude and Longitude Vectors
% Define latitude vector in descending order (from north to south)
pop_lat_vec_full = linspace(pop_R.LatitudeLimits(2), pop_R.LatitudeLimits(1), size(pop_data_full,1))';
% Define longitude vector in ascending order (from west to east)
pop_lon_vec_full = linspace(pop_R.LongitudeLimits(1), pop_R.LongitudeLimits(2), size(pop_data_full,2));

%% 5. Apply Offsets to Align Grid Cells
% Define a small offset based on the step size for latitude and longitude
lat_step = abs(pop_lat_vec_full(2) - pop_lat_vec_full(1));
lon_step = abs(pop_lon_vec_full(2) - pop_lon_vec_full(1));
lat_offset = lat_step * 5; % Part of a cell size for latitude
lon_offset = lon_step / 2.5; % Part of a cell size for longitude

% Apply the offset to shift the cells appropriately
% Since latitude is descending, subtracting the offset shifts it downwards
pop_lat_vec_full = pop_lat_vec_full - lat_offset;
% Adding the offset to longitude shifts it to the right
pop_lon_vec_full = pop_lon_vec_full + lon_offset;

%% 6. Extract Subset of Population Data Corresponding to Region of Interest
% Find indices for latitude
lat_indices = find(pop_lat_vec_full >= lat_min & pop_lat_vec_full <= lat_max);
% Find indices for longitude
lon_indices = find(pop_lon_vec_full >= lon_min & pop_lon_vec_full <= lon_max);

% Extract the subset of data
pop_data = pop_data_full(lat_indices, lon_indices);
pop_lat_vec = pop_lat_vec_full(lat_indices);
pop_lon_vec = pop_lon_vec_full(lon_indices);

% Create latitude and longitude grids
[pop_lon_grid, pop_lat_grid] = meshgrid(pop_lon_vec, pop_lat_vec);

%% 7. Correct Data Orientation Based on Spatial Referencing
% Flip vertically if necessary
if strcmp(pop_R.ColumnsStartFrom, 'north')
    pop_data = flipud(pop_data);
    disp('Flipped population data vertically to align latitude from south to north.');
end

%% 8. Visualize Population Data with Proper Scaling
figure;
imagesc(pop_lon_vec, pop_lat_vec, log10(pop_data + 1)); % Log scale for better visibility
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Population Data (10³) - Region of Interest');
% Set fixed aspect ratio
axis image;

%% 10. Visualize Non-Zero Areas
figure;
imagesc(pop_lon_vec, pop_lat_vec, pop_data > 0);
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Areas with Non-Zero Population - Region of Interest');
% Set fixed aspect ratio
axis image;

%% 11. Aggregate Population Data onto CFD Grid

% Calculate step sizes based on grid centers
step_lon = (lon_max - lon_min) / GridX;
step_lat = (lat_max - lat_min) / GridY;

% Define CFD grid edges based on cell centers
lon_edges = linspace(lon_min - step_lon/2, lon_max + step_lon/2, GridX + 1);
lat_edges = linspace(lat_min - step_lat/2, lat_max + step_lat/2, GridY + 1);

% Flatten the population data grids
pop_lat_flat = pop_lat_grid(:);
pop_lon_flat = pop_lon_grid(:);
pop_data_flat = pop_data(:);

% Assign Each Data Point to a CFD Grid Cell using bin indices
[~, ~, lon_bin_idx] = histcounts(pop_lon_flat, lon_edges);
[~, ~, lat_bin_idx] = histcounts(pop_lat_flat, lat_edges);

% Remove data points that fall outside the CFD grid (bin index 0)
valid_idx = lat_bin_idx > 0 & lon_bin_idx > 0;
lat_bin_idx = lat_bin_idx(valid_idx);
lon_bin_idx = lon_bin_idx(valid_idx);
pop_data_flat = pop_data_flat(valid_idx);

% Verify that all valid data points are within grid bounds
assert(all(lat_bin_idx >= 1 & lat_bin_idx <= GridY), 'Latitude bin indices out of range.');
assert(all(lon_bin_idx >= 1 & lon_bin_idx <= GridX), 'Longitude bin indices out of range.');

% Aggregate population data into CFD grid cells
population_on_CFD_grid = accumarray([lat_bin_idx, lon_bin_idx], pop_data_flat, [GridY, GridX], @sum, 0);

% Convert the result to integer type for memory efficiency
population_on_CFD_grid = int32(round(population_on_CFD_grid));

% Verify that the total population remains the same after aggregation
original_total_pop = sum(pop_data_flat);
aggregated_total_pop = sum(population_on_CFD_grid(:));

tolerance = 1000; % Allowable difference between original and aggregated total
assert(abs(original_total_pop - aggregated_total_pop) <= tolerance, ...
       ['Total population mismatch: Original = ' num2str(original_total_pop) ...
        ', Aggregated = ' num2str(aggregated_total_pop)]);

disp(['Total Population (Original): ', num2str(original_total_pop)]);
disp(['Total Population (Aggregated): ', num2str(aggregated_total_pop)]);

% Generate CFD grid coordinates for visualization
[X_CFD_edges, Y_CFD_edges] = meshgrid(lon_edges, lat_edges);

% Generate center points of CFD grid cells for plotting
X_CFD_centers = (X_CFD_edges(1:end-1,1:end-1) + X_CFD_edges(1:end-1,2:end)) / 2;
Y_CFD_centers = (Y_CFD_edges(1:end-1,1:end-1) + Y_CFD_edges(2:end,1:end-1)) / 2;

% Display statistics of aggregated data
min_agg = min(population_on_CFD_grid(:));
max_agg = max(population_on_CFD_grid(:));
mean_agg = mean(population_on_CFD_grid(:));
disp(['Aggregated Population Data - Min: ', num2str(min_agg), ...
      ', Max: ', num2str(max_agg), ', Mean: ', num2str(mean_agg)]);

population_on_CFD_grid = double(population_on_CFD_grid);

%% 12. Visualize Aggregated Population Data
figure;
imagesc(lon_edges, lat_edges, log10(population_on_CFD_grid + 1)); % Log scale
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Aggregated Population on CFD Grid (10³)');
% Set fixed aspect ratio
% axis image;%% Define the Specific Grid Area to Check
y_range = 64:65;
x_range = 95:96;

%% 1. Get Limits for the Area in Latitude and Longitude
% Calculate the latitude and longitude limits of the specified area
lat_limits_check = [lat_edges(y_range(1)), lat_edges(y_range(end) + 1)];
lon_limits_check = [lon_edges(x_range(1)), lon_edges(x_range(end) + 1)];

% Display the latitude and longitude limits
disp('Latitude Limits of the Area to Check:');
disp(lat_limits_check);
disp('Longitude Limits of the Area to Check:');
disp(lon_limits_check);

%% 2. Find Indices in the Original Population Data Based on Latitude and Longitude
% Find indices in the latitude vector for the limits
lat_indices_check = find(pop_lat_vec_full >= lat_limits_check(1) & pop_lat_vec_full <= lat_limits_check(2));
% Find indices in the longitude vector for the limits
lon_indices_check = find(pop_lon_vec_full >= lon_limits_check(1) & pop_lon_vec_full <= lon_limits_check(2));

% Ensure that indices are within the bounds of the original data grid
lat_indices_check = lat_indices_check(lat_indices_check > 0 & lat_indices_check <= size(pop_data_full, 1));
lon_indices_check = lon_indices_check(lon_indices_check > 0 & lon_indices_check <= size(pop_data_full, 2));

% Extract the original population data subset based on latitude and longitude
original_pop_subset = pop_data_full(lat_indices_check, lon_indices_check);
original_pop_sum = round(sum(original_pop_subset(:)));

%% 3. Extract Population Data from the Aggregated CFD Grid for the Specified Area
% Extract aggregated population data for the specified area
aggregated_pop_subset = population_on_CFD_grid(y_range, x_range);
aggregated_pop_sum = round(sum(aggregated_pop_subset(:)));

% Display Results
disp('Original Population Sum for Specified Area:');
disp(original_pop_sum);
disp('Aggregated Population Sum for Specified Area:');
disp(aggregated_pop_sum);

%% 4. Zoomed Population Data from the Aggregated CFD Grid
% Define the x and y ranges for the CFD grid to zoom into
% Define the x and y ranges for the CFD grid to zoom into
y_range = 60:70; % X grid indices to display
x_range = 90:100; % Y grid indices to display

% Define the zoomed region limits
zoom_lat_min = lat_edges(y_range(1));
zoom_lat_max = lat_edges(y_range(end));
zoom_lon_min = lon_edges(x_range(1));
zoom_lon_max = lon_edges(x_range(end));


% Original Population data
figure;

% Set up map projection for the zoomed region
axesm('mercator', 'MapLatLimit', [zoom_lat_min, zoom_lat_max], 'MapLonLimit', [zoom_lon_min, zoom_lon_max]);
axis off;

% Find an OpenStreetMap WMS layer to use for the background map
osmLayer = wmsfind('osm', 'SearchField', 'serverurl');  
osmLayer = refine(osmLayer, 'osm');
layer = osmLayer(1);

% Retrieve the map data from the WMS server
imageWidth = 2000;  
imageHeight = 2000;

% Retrieve zoomed map data from the WMS server
[A_zoom, R_zoom] = wmsread(layer, 'Latlim', [zoom_lat_min, zoom_lat_max], 'Lonlim', [zoom_lon_min, zoom_lon_max], ...
    'ImageHeight', imageHeight, 'ImageWidth', imageWidth);

hMap1 = imshow(A_zoom, 'XData', [zoom_lon_min zoom_lon_max], 'YData', [zoom_lat_min zoom_lat_max]);  
set(hMap1, 'AlphaData', 0.5);

hold on;

% Overlay the population data and grid
cont1 = imagesc(pop_lon_vec, pop_lat_vec, log10(round(pop_data))); % Adjust AlphaData for transparency
uistack(cont1, 'bottom');
colorbar;
% Set fixed aspect ratio
axis image;
% Set axis limits to zoom in on the region of interest
xlim([zoom_lon_min, zoom_lon_max]);
ylim([zoom_lat_min, zoom_lat_max]);
title('Zoomed Overlay of Population Data and CFD Grid Boundaries with Map Background');
hold off;

%% CFD Aggregated

figure;

% Set up map projection for the zoomed region
axesm('mercator', 'MapLatLimit', [zoom_lat_min, zoom_lat_max], 'MapLonLimit', [zoom_lon_min, zoom_lon_max]);
axis off;

% Find an OpenStreetMap WMS layer to use for the background map
osmLayer = wmsfind('osm', 'SearchField', 'serverurl');  
osmLayer = refine(osmLayer, 'osm');
layer = osmLayer(1);

% Retrieve the map data from the WMS server
imageWidth = 2000;  
imageHeight = 2000;

% Retrieve zoomed map data from the WMS server
[A_zoom, R_zoom] = wmsread(layer, 'Latlim', [zoom_lat_min, zoom_lat_max], 'Lonlim', [zoom_lon_min, zoom_lon_max], ...
    'ImageHeight', imageHeight, 'ImageWidth', imageWidth);

hMap1 = imshow(A_zoom, 'XData', [zoom_lon_min zoom_lon_max], 'YData', [zoom_lat_min zoom_lat_max]);  
set(hMap1, 'AlphaData', 0.5);

hold on;
% lon_edges, lat_edges, log10(population_on_CFD_grid + 1)
% Overlay the population data and grid
cont1 = imagesc(lon_edges-step_lon/2, lat_edges+step_lat/2, log10(round(population_on_CFD_grid+1))); % Offseted so the plot is in the center of the cell
uistack(cont1, 'bottom');
colorbar;
% Set fixed aspect ratio
axis image;
% Set axis limits to zoom in on the region of interest
xlim([zoom_lon_min, zoom_lon_max]);
ylim([zoom_lat_min, zoom_lat_max]);
title('Zoomed Overlay of Population Data and CFD Grid Boundaries with Map Background');
hold off;

%% 13. Save the Aggregated Population Data
save('population_on_CFD_grid.mat', 'population_on_CFD_grid', ...
     'X_CFD_centers', 'Y_CFD_centers');

disp('Population data has been aggregated onto the CFD grid and saved.');