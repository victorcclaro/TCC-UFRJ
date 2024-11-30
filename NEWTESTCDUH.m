clear all; close all; clc;

% -----------------------------
% Load Necessary Parameters
% -----------------------------

height = 100;  % Change this to the desired height
num_years = 5;

% Load the CFD results including snapshots
cfd_data_filename = sprintf('cfd_results_%dm_combined_%dyears.mat', height, num_years);
load(cfd_data_filename, 'Nt', 'deposition_total', 'average_concentration_kgm3', 'dC_dM', ...
    'delta_x', 'delta_y', 'delta_t', 'lat_min', 'lat_max', ...
    'lon_min', 'lon_max', 'GridX', 'GridY');

% Load interpolated population data for calculations and visualization
load('4newpopulation_on_CFD_grid.mat', 'population_on_CFD_grid', 'X_CFD_centers', 'Y_CFD_centers');



%% -----------------------------
% Health Impact Calculations
% -----------------------------

% Step 1: Data Preparation

% Convert average concentration from kg/m³ to µg/m³
average_concentration_mugm3 = average_concentration_kgm3 * 1e9;  % µg/m³

% Threshold for impacted cells (set as needed)
threshold_concentration = 1;  % µg/m³

% Identify impacted cells based on concentration threshold
impacted_cells = average_concentration_mugm3 >= threshold_concentration;

% Step 2: Define Health Impact Parameters

% Health Effects:
% 1 - Chronic Mortality
% 2 - Acute Mortality
% 3 - Respiratory Morbidity
% 4 - Cardiovascular Morbidity

% Relative Risks (RR) per µg/m³ 
RR = [1.0043, 1.0006, 1.00114, 1.0005];

% Incidence Rates (Finc,e) per year per person 
Finc = [6.76e-3, 6.76e-3, 3.08e-3, 5.28e-3];

% Years of Life Lost (YLLe) 
YLLe = [10, 0.25, 0, 0];  % [Chronic, Acute, Respiratory, Cardiovascular]

% Duration of Health Effect (De) in years 
De = [0, 0, 0.04, 0.04];  % [Chronic, Acute, Respiratory, Cardiovascular]

% Severity Factor (Se) 
Se = [0, 0, 0.64, 0.71];  % [Chronic, Acute, Respiratory, Cardiovascular]

% Population per grid cell
Ni = population_on_CFD_grid';  % Transposed to match grid orientation

% Concentration per grid cell
Ck_i = average_concentration_mugm3;  % 2D matrix, same size

% Step 3: Initialize DALYs Matrices

% Initialize DALYs matrices (same size as average_concentration_mugm3)
DALYs_chronic_mortality = zeros(size(Ck_i));
DALYs_acute_mortality = zeros(size(Ck_i));
DALYs_respiratory_morbidity = zeros(size(Ck_i));
DALYs_cardiovascular_morbidity = zeros(size(Ck_i));

% Step 4: Calculate ABe and DALYs for Each Health Effect

% 1. Chronic Mortality
ABe_chronic = ((RR(1) - 1) .* Ck_i .* Finc(1)) ./ ((RR(1) - 1) .* Ck_i + 1);
Cases_chronic_mortality = ABe_chronic .* Ni;
DALYs_chronic_mortality = Cases_chronic_mortality .* YLLe(1);

% 2. Acute Mortality
ABe_acute = ((RR(2) - 1) .* Ck_i .* Finc(2)) ./ ((RR(2) - 1) .* Ck_i + 1);
Cases_acute_mortality = ABe_acute .* Ni;
DALYs_acute_mortality = Cases_acute_mortality .* YLLe(2);

% 3. Respiratory Morbidity
ABe_respiratory = ((RR(3) - 1) .* Ck_i .* Finc(3)) ./ ((RR(3) - 1) .* Ck_i + 1);
Cases_respiratory_morbidity = ABe_respiratory .* Ni;
DALYs_respiratory_morbidity = Cases_respiratory_morbidity .* (De(3) * Se(3));

% 4. Cardiovascular Morbidity
ABe_cardiovascular = ((RR(4) - 1) .* Ck_i .* Finc(4)) ./ ((RR(4) - 1) .* Ck_i + 1);
Cases_cardiovascular_morbidity = ABe_cardiovascular .* Ni;
DALYs_cardiovascular_morbidity = Cases_cardiovascular_morbidity .* (De(4) * Se(4));

% Step 5: Sum DALYs per Grid Cell

Cases_total_morbidity =  Cases_respiratory_morbidity + Cases_cardiovascular_morbidity;
Cases_total_mortality = Cases_chronic_mortality + Cases_acute_mortality;
DALYs_total = DALYs_chronic_mortality + DALYs_acute_mortality + ...
             DALYs_respiratory_morbidity + DALYs_cardiovascular_morbidity;

% Step 6: Apply Concentration Threshold

% Set Cases to zero for grid cells below the concentration threshold
Cases_chronic_mortality(~impacted_cells) = 0;
Cases_acute_mortality(~impacted_cells) = 0;
Cases_respiratory_morbidity(~impacted_cells) = 0;
Cases_cardiovascular_morbidity(~impacted_cells) = 0;
Cases_total_morbidity(~impacted_cells) = 0;
Cases_total_mortality(~impacted_cells) = 0;

DALYs_chronic_mortality(~impacted_cells) = 0;
DALYs_acute_mortality(~impacted_cells) = 0;
DALYs_respiratory_morbidity(~impacted_cells) = 0;
DALYs_cardiovascular_morbidity(~impacted_cells) = 0;
DALYs_total(~impacted_cells) = 0;

% Step 7: Calculate Total Health Impacts Over Impacted Cells

% Initialize VR Cases and DALYs
VR_Cases_chronic_mortality_impacted = 0;
VR_Cases_acute_mortality_impacted = 0;
VR_Cases_respiratory_morbidity_impacted = 0;
VR_Cases_cardiovascular_morbidity_impacted = 0;

VR_DALYs_chronic_mortality_impacted = 0;
VR_DALYs_acute_mortality_impacted = 0;
VR_DALYs_respiratory_morbidity_impacted = 0;
VR_DALYs_cardiovascular_morbidity_impacted = 0;

population_VR = 0;

% Define VR-specific grid indices (corrected to include distinct cells)
vr_indices = [
    95, 64;
    96, 64;
    95, 65;
    96, 65
];

% Sum Cases and DALYs for each VR grid cell
for i = 1:size(vr_indices,1)
    row = vr_indices(i,1);
    col = vr_indices(i,2);
    
    VR_Cases_chronic_mortality_impacted = VR_Cases_chronic_mortality_impacted + Cases_chronic_mortality(row, col);
    VR_Cases_acute_mortality_impacted = VR_Cases_acute_mortality_impacted + Cases_acute_mortality(row, col);
    VR_Cases_respiratory_morbidity_impacted = VR_Cases_respiratory_morbidity_impacted + Cases_respiratory_morbidity(row, col);
    VR_Cases_cardiovascular_morbidity_impacted = VR_Cases_cardiovascular_morbidity_impacted + Cases_cardiovascular_morbidity(row, col);
    
    VR_DALYs_chronic_mortality_impacted = VR_DALYs_chronic_mortality_impacted + DALYs_chronic_mortality(row, col);
    VR_DALYs_acute_mortality_impacted = VR_DALYs_acute_mortality_impacted + DALYs_acute_mortality(row, col);
    VR_DALYs_respiratory_morbidity_impacted = VR_DALYs_respiratory_morbidity_impacted + DALYs_respiratory_morbidity(row, col);
    VR_DALYs_cardiovascular_morbidity_impacted = VR_DALYs_cardiovascular_morbidity_impacted + DALYs_cardiovascular_morbidity(row, col);

    population_VR = population_VR + population_on_CFD_grid(col, row);
end

% Calculate VR Health Impacts
VR_Cases_chronic_mortality_impacted = round(VR_Cases_chronic_mortality_impacted, 2);
VR_Cases_acute_mortality_impacted = round(VR_Cases_acute_mortality_impacted, 2);
VR_Cases_respiratory_morbidity_impacted = round(VR_Cases_respiratory_morbidity_impacted, 2);
VR_Cases_cardiovascular_morbidity_impacted = round(VR_Cases_cardiovascular_morbidity_impacted, 2);
VR_Cases_mortality_impacted = round(VR_Cases_chronic_mortality_impacted + VR_Cases_acute_mortality_impacted, 2);
VR_Cases_morbidity_impacted = round(VR_Cases_respiratory_morbidity_impacted + VR_Cases_cardiovascular_morbidity_impacted, 2);
VR_Total_Cases_Impacted = round(VR_Cases_mortality_impacted + VR_Cases_morbidity_impacted, 2);

% VR DALYs 
VR_DALYs_chronic_mortality_impacted = round(VR_DALYs_chronic_mortality_impacted, 2);
VR_DALYs_acute_mortality_impacted = round(VR_DALYs_acute_mortality_impacted, 2);
VR_DALYs_respiratory_morbidity_impacted = round(VR_DALYs_respiratory_morbidity_impacted, 2);
VR_DALYs_cardiovascular_morbidity_impacted = round(VR_DALYs_cardiovascular_morbidity_impacted, 2);
VR_DALYs_mortality_impacted = round(VR_DALYs_chronic_mortality_impacted + VR_DALYs_acute_mortality_impacted, 2);
VR_DALYs_morbidity_impacted = round(VR_DALYs_respiratory_morbidity_impacted + VR_DALYs_cardiovascular_morbidity_impacted, 2);
VR_Total_DALYs_Impacted = round(VR_DALYs_mortality_impacted + VR_DALYs_morbidity_impacted, 2);
VR_Total_DALYs_Impacted_perPerson_perLifetime = VR_Total_DALYs_Impacted * 80 / 100000;

% Calculate Total Health Impacts Across All Impacted Cells
total_Cases_chronic_mortality_impacted = round(sum(Cases_chronic_mortality(:)), 2);
total_Cases_acute_mortality_impacted = round(sum(Cases_acute_mortality(:)), 2);
total_Cases_respiratory_morbidity_impacted = round(sum(Cases_respiratory_morbidity(:)), 2);
total_Cases_cardiovascular_morbidity_impacted = round(sum(Cases_cardiovascular_morbidity(:)), 2);
total_Cases_mortality_impacted = round(total_Cases_chronic_mortality_impacted + total_Cases_acute_mortality_impacted, 2);
total_Cases_morbidity_impacted = round(total_Cases_respiratory_morbidity_impacted + total_Cases_cardiovascular_morbidity_impacted, 2);
total_Cases_impacted = round(total_Cases_morbidity_impacted(:) + total_Cases_mortality_impacted(:),2);

% Total DALYs Across All Impacted Cells
total_DALYs_impacted = round(sum(DALYs_total(:)), 2);
total_DALYs_chronic_mortality_impacted = round(sum(DALYs_chronic_mortality(:)), 2);
total_DALYs_acute_mortality_impacted = round(sum(DALYs_acute_mortality(:)), 2);
total_DALYs_respiratory_morbidity_impacted = round(sum(DALYs_respiratory_morbidity(:)), 2);
total_DALYs_cardiovascular_morbidity_impacted = round(sum(DALYs_cardiovascular_morbidity(:)), 2);
total_DALYs_mortality_impacted = round(total_DALYs_chronic_mortality_impacted + total_DALYs_acute_mortality_impacted, 2);
total_DALYs_morbidity_impacted = round(total_DALYs_respiratory_morbidity_impacted + total_DALYs_cardiovascular_morbidity_impacted, 2);

% Calculate Percentages for VR Cases and DALYs Relative to Total Impacted Region

percent_Cases_respiratory_morbidity_VR = (VR_Cases_respiratory_morbidity_impacted / total_Cases_respiratory_morbidity_impacted) * 100;
percent_Cases_cardiovascular_morbidity_VR = (VR_Cases_cardiovascular_morbidity_impacted / total_Cases_cardiovascular_morbidity_impacted) * 100;
percent_Cases_morbidity_VR = (VR_Cases_morbidity_impacted / total_Cases_morbidity_impacted) * 100;
percent_Cases_chronic_mortality_VR = (VR_Cases_chronic_mortality_impacted / total_Cases_chronic_mortality_impacted) * 100;
percent_Cases_acute_mortality_VR = (VR_Cases_acute_mortality_impacted / total_Cases_acute_mortality_impacted) * 100;
percent_Cases_mortality_VR = (VR_Cases_mortality_impacted / total_Cases_mortality_impacted) * 100;

percent_DALYs_respiratory_morbidity_VR = (VR_DALYs_respiratory_morbidity_impacted / total_DALYs_respiratory_morbidity_impacted) * 100;
percent_DALYs_cardiovascular_morbidity_VR = (VR_DALYs_cardiovascular_morbidity_impacted / total_DALYs_cardiovascular_morbidity_impacted) * 100;
percent_DALYs_morbidity_VR = (VR_DALYs_morbidity_impacted / total_DALYs_morbidity_impacted) * 100;
percent_DALYs_chronic_mortality_VR = (VR_DALYs_chronic_mortality_impacted / total_DALYs_chronic_mortality_impacted) * 100;
percent_DALYs_acute_mortality_VR = (VR_DALYs_acute_mortality_impacted / total_DALYs_acute_mortality_impacted) * 100;
percent_DALYs_mortality_VR = (VR_DALYs_mortality_impacted / total_DALYs_mortality_impacted) * 100;

% -----------------------------
% Step 8: Display Results
% -----------------------------

fprintf('* Health Impact Results *\n\n');

Relative_factor = 100000/population_VR;

fprintf('* VR Results per 100,000 people*\n');
fprintf('--- VR Cases Results ---\n');
fprintf('Respiratory Morbidity per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_Cases_respiratory_morbidity_impacted*Relative_factor, percent_Cases_respiratory_morbidity_VR);
fprintf('Cardiovascular Morbidity per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_Cases_cardiovascular_morbidity_impacted*Relative_factor, percent_Cases_cardiovascular_morbidity_VR);
fprintf('* VR Total Cases from Morbidity per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_Cases_morbidity_impacted*Relative_factor, percent_Cases_morbidity_VR);
fprintf('Chronic Mortality per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_Cases_chronic_mortality_impacted*Relative_factor, percent_Cases_chronic_mortality_VR);
fprintf('Acute Mortality per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_Cases_acute_mortality_impacted*Relative_factor, percent_Cases_acute_mortality_VR);
fprintf('* VR Total Cases from Mortality per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_Cases_mortality_impacted*Relative_factor, percent_Cases_mortality_VR);
fprintf('** VR Total Cases Impacted per 100,000 people: %.2f\n\n', VR_Total_Cases_Impacted*Relative_factor);

fprintf('--- VR DALYs Results ---\n');
fprintf('DALYs from Respiratory Morbidity: %.2f (%.2f%% of impacted region)\n', ...
    VR_DALYs_respiratory_morbidity_impacted*Relative_factor, percent_DALYs_respiratory_morbidity_VR);
fprintf('DALYs from Cardiovascular Morbidity per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_DALYs_cardiovascular_morbidity_impacted*Relative_factor, percent_DALYs_cardiovascular_morbidity_VR);
fprintf('* VR Total DALYs from Morbidity per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_DALYs_morbidity_impacted*Relative_factor, percent_DALYs_morbidity_VR);
fprintf('DALYs from Chronic Mortality per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_DALYs_chronic_mortality_impacted*Relative_factor, percent_DALYs_chronic_mortality_VR);
fprintf('DALYs from Acute Mortality per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_DALYs_acute_mortality_impacted*Relative_factor, percent_DALYs_acute_mortality_VR);
fprintf('* VR Total DALYs from Mortality per 100,000 people: %.2f (%.2f%% of impacted region)\n', ...
    VR_DALYs_mortality_impacted*Relative_factor, percent_DALYs_mortality_VR);
fprintf('** VR Total DALYs Impacted per 100,000 people: %.2f\n\n', VR_Total_DALYs_Impacted*Relative_factor);
fprintf('** VR Total DALYs Impacted per person per lifetime (80y): %.2f\n\n', VR_Total_DALYs_Impacted_perPerson_perLifetime*Relative_factor);

fprintf('* All Impacted Area Results *\n');
fprintf('--- Cases Results ---\n');
fprintf('Respiratory Morbidity: %.2f\n', total_Cases_respiratory_morbidity_impacted);
fprintf('Cardiovascular Morbidity: %.2f\n', total_Cases_cardiovascular_morbidity_impacted);
fprintf('* Total Cases from Morbidity: %.2f\n', total_Cases_morbidity_impacted);
fprintf('Chronic Mortality: %.2f\n', total_Cases_chronic_mortality_impacted);
fprintf('Acute Mortality: %.2f\n', total_Cases_acute_mortality_impacted);
fprintf('* Total Cases from Mortality: %.2f\n', total_Cases_mortality_impacted);
fprintf('** Total Cases Impacted: %.2f\n\n', total_Cases_impacted);

fprintf('--- DALYs Results ---\n');
fprintf('Respiratory Morbidity: %.2f\n', total_DALYs_respiratory_morbidity_impacted);
fprintf('Cardiovascular Morbidity: %.2f\n', total_DALYs_cardiovascular_morbidity_impacted);
fprintf('* Total DALYs from Morbidity: %.2f\n', total_DALYs_morbidity_impacted);
fprintf('Chronic Mortality: %.2f\n', total_DALYs_chronic_mortality_impacted);
fprintf('Acute Mortality: %.2f\n', total_DALYs_acute_mortality_impacted);
fprintf('* Total DALYs from Mortality: %.2f\n', total_DALYs_mortality_impacted);
fprintf('** Total DALYs Impacted: %.2f\n', total_DALYs_impacted);

%% -----------------------------
% Step 8: Calculate deposition and impacted population
% -----------------------------

% Total simulation time in seconds
total_time_sec = Nt * delta_t;  

% Scale the total deposition to per year
deposition_total_annual = deposition_total * (365 * 24 * 3600 / total_time_sec);  % kg/m²/year

% Calculate the area of a grid cell
cell_area = delta_x * delta_y;  % m²

% Calculate total deposition per grid cell
deposition_total_cell_annual = deposition_total_annual * cell_area;  % kg/year per cell

% Calculate the total deposition over impacted cells
total_deposition_impacted_kg = round(sum(deposition_total_cell_annual(impacted_cells)));
fprintf('Total deposition in impacted area (kg): %d\n', total_deposition_impacted_kg);

% Calculate the total impacted population
total_population_impacted = round(sum(population_on_CFD_grid(impacted_cells')));
fprintf('Impacted population: %d\n', total_population_impacted);
fprintf('Population on VR: %d\n', round(population_VR));

%% -----------------------------
% Adjust max and min latitudes and longitudes to impacted area
% -----------------------------

[LonMesh, LatMesh] = meshgrid(linspace(lon_min, lon_max, GridX), linspace(lat_min, lat_max, GridY));

% Filter latitude and longitude for impacted cells based on logical indexing
latitudes_impacted = LatMesh(impacted_cells');
longitudes_impacted = LonMesh(impacted_cells');

% Calculate min and max latitude and longitude for impacted cells
lat_min_impacted = min(latitudes_impacted)-0.5;
lat_max_impacted = max(latitudes_impacted)+0.1;
lon_min_impacted = min(longitudes_impacted)-0.5;
lon_max_impacted = max(longitudes_impacted)+0.5;

% Display the results
fprintf('Latitude range for impacted area: [%f, %f]\n', lat_min_impacted, lat_max_impacted);
fprintf('Longitude range for impacted area: [%f, %f]\n', lon_min_impacted, lon_max_impacted);

%% -----------------------------
% Visualization of Health Results
% -----------------------------

% Retrieve the Map Image from WMS Service
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
scale_length_km = 50;  % 50 km scale bar

% Calculate the scale bar length in degrees of longitude
scale_bar_length_deg = scale_length_km * degrees_per_km_lon;

% Define offsets from the map edges (as a percentage of map dimensions)
offset_x_pct = 0.05;  % 5% offset from the right edge
offset_y_pct = 0.05;  % 5% offset from the bottom edge

% Calculate the offsets in degrees
offset_x_deg = (lon_max_impacted - lon_min_impacted) * offset_x_pct;
offset_y_deg = (lat_max_impacted - lat_min_impacted) * offset_y_pct;

% Define the starting point of the scale bar
x_start = lon_max_impacted - offset_x_deg - scale_bar_length_deg;
y_start = lat_min_impacted + offset_y_deg;

% Define the end point of the scale bar
x_end = x_start + scale_bar_length_deg;
y_end = y_start;

% Generate Meshgrid for Contour Plots
[X, Y] = meshgrid(linspace(lon_min, lon_max, GridX), linspace(lat_min, lat_max, GridY));

% Create a new figure for visualization
figure (1);

% Subplot 1: Average PM10 Concentration (µg/m³)
subplot(1, 2, 1);
hMap1 = imshow(A, 'XData', [lon_min lon_max], 'YData', [lat_min lat_max]);  
set(hMap1, 'AlphaData', 0.5);
axis on;        
axis equal;     
hold on;  

% Plot the scale bar
plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', 2);
plot([x_start, x_start], [y_start - offset_y_deg * 0.1, y_start + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
plot([x_end, x_end], [y_end - offset_y_deg * 0.1, y_end + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
text((x_start + x_end)/2, y_start - 0.5*offset_y_deg, ...
    [num2str(scale_length_km) ' km'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'FontWeight', 'bold');

% Create the contour plot for PM10 concentration
[c1, cont1] = contourf(X, Y, average_concentration_mugm3', 'LineColor', 'none');
set(gca, 'XLim', [lon_min_impacted lon_max_impacted], 'YLim', [lat_min_impacted lat_max_impacted]);
colorbar;
uistack(cont1, 'bottom');

% Titles and labels
title('Average PM_{10} Concentration (\mug/m^3)');
xlabel('Longitude');
ylabel('Latitude');

% Highlight impacted cells
contour(X, Y, impacted_cells', [1 1], 'LineColor', 'r', 'LineWidth', 2);

% Subplot 2: Total DALYs per Cell
subplot(1, 2, 2);
hMap2 = imshow(A, 'XData', [lon_min lon_max], 'YData', [lat_min lat_max]);  
set(hMap2, 'AlphaData', 0.5);
axis on;        
axis equal;     
hold on;  

% Plot the scale bar
plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', 2);
plot([x_start, x_start], [y_start - offset_y_deg * 0.1, y_start + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
plot([x_end, x_end], [y_end - offset_y_deg * 0.1, y_end + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
text((x_start + x_end)/2, y_start - 0.5*offset_y_deg, ...
    [num2str(scale_length_km) ' km'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'FontWeight', 'bold');

% Create the contour plot for DALYs
[c2, cont2] = contourf(X, Y, DALYs_total', 'LineColor', 'none');
set(gca, 'XLim', [lon_min_impacted lon_max_impacted], 'YLim', [lat_min_impacted lat_max_impacted]);
colorbar;
uistack(cont2, 'bottom');

% Titles and labels
title('Total DALYs per Year');
xlabel('Longitude');
ylabel('Latitude');


% Display total DALYs
text(lon_max_impacted - offset_x_deg, lat_max_impacted - 1.5*offset_y_deg, ...
    sprintf('Total DALYs:\n%.0f years', total_DALYs_impacted), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'FontWeight', 'bold', 'BackgroundColor', 'w');

% Highlight impacted cells
contour(X, Y, impacted_cells', [1 1], 'LineColor', 'r', 'LineWidth', 2);

% Create a new figure for visualization
figure (2);

% Subplot 3: Total Mortality Cases per Cell
subplot(1, 2, 1);
hMap3 = imshow(A, 'XData', [lon_min lon_max], 'YData', [lat_min lat_max]);  
set(hMap3, 'AlphaData', 0.5);
axis on;        
axis equal;     
hold on;  

% Plot the scale bar
plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', 2);
plot([x_start, x_start], [y_start - offset_y_deg * 0.1, y_start + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
plot([x_end, x_end], [y_end - offset_y_deg * 0.1, y_end + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
text((x_start + x_end)/2, y_start - 0.5*offset_y_deg, ...
    [num2str(scale_length_km) ' km'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'FontWeight', 'bold');

% Create the contour plot for mortality cases
[c3, cont3] = contourf(X, Y, Cases_total_mortality', 'LineColor', 'none');
set(gca, 'XLim', [lon_min_impacted lon_max_impacted], 'YLim', [lat_min_impacted lat_max_impacted]);
colorbar;
uistack(cont3, 'bottom');

% Titles and labels
title('Annual Mortality Cases per Year');
xlabel('Longitude');
ylabel('Latitude');

% Display total mortality cases
text(lon_max_impacted - offset_x_deg, lat_max_impacted - 1.5*offset_y_deg, ...
    sprintf('Total Mortality:\n%.0f cases', total_Cases_mortality_impacted), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'FontWeight', 'bold', 'BackgroundColor', 'w');

% Highlight impacted cells
contour(X, Y, impacted_cells', [1 1], 'LineColor', 'r', 'LineWidth', 2);

% Subplot 4: Total Morbidity Cases per Cell (Hospitalizations)
subplot(1, 2, 2);
hMap4 = imshow(A, 'XData', [lon_min lon_max], 'YData', [lat_min lat_max]);  
set(hMap4, 'AlphaData', 0.5);
axis on;        
axis equal;     
hold on;  

% Plot the scale bar
plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', 2);
plot([x_start, x_start], [y_start - offset_y_deg * 0.1, y_start + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
plot([x_end, x_end], [y_end - offset_y_deg * 0.1, y_end + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
text((x_start + x_end)/2, y_start - offset_y_deg, ...
    [num2str(scale_length_km) ' km'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'FontWeight', 'bold');

% Create the contour plot for morbidity cases
[c4, cont4] = contourf(X, Y, Cases_total_morbidity', 'LineColor', 'none');
set(gca, 'XLim', [lon_min_impacted lon_max_impacted], 'YLim', [lat_min_impacted lat_max_impacted]);
colorbar;
uistack(cont4, 'bottom');

% Titles and labels
title('Annual Morbidity Cases per Year');
xlabel('Longitude');
ylabel('Latitude');

% Display total morbidity cases
text(lon_max_impacted - offset_x_deg, lat_max_impacted - 1.5*offset_y_deg, ...
    sprintf('Total Morbidity:\n%.0f cases', total_Cases_morbidity_impacted), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'FontWeight', 'bold', 'BackgroundColor', 'w');

% Highlight impacted cells
contour(X, Y, impacted_cells', [1 1], 'LineColor', 'r', 'LineWidth', 2);

%% -----------------------------
% Visualization of Deposition and Population
% -----------------------------

figure (3);

% Subplot 5: Total Annual Deposition per Cell
% Display the map with geographic referencing
subplot(1, 2, 1);
hMap1 = imshow(A, 'XData', [lon_min lon_max], 'YData', [lat_min lat_max]);  
set(hMap1, 'AlphaData', 0.5);

% Ensure that the aspect ratio is correct
axis on;        % Turn on axis
axis equal;     % Keep the geographic aspect ratio

hold on;  % Hold on to overlay the contour plot

% Plot the scale bar
plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', 2);
plot([x_start, x_start], [y_start - offset_y_deg * 0.1, y_start + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
plot([x_end, x_end], [y_end - offset_y_deg * 0.1, y_end + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
text(x_start + scale_bar_length_deg/2, y_start - 0.5*offset_y_deg, ...
    [num2str(scale_length_km) ' km'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'FontWeight', 'bold');

% Generate Meshgrid for Contour Plot
[X, Y] = meshgrid(linspace(lon_min, lon_max, GridX), linspace(lat_min, lat_max, GridY));

% Create the contour plot
[c1, cont1] = contourf(X, Y, deposition_total_cell_annual', 'LineColor', 'none');
set(gca, 'XLim', [lon_min_impacted lon_max_impacted], 'YLim', [lat_min_impacted lat_max_impacted]);
colorbar;
uistack(cont1, 'bottom');

% Titles and labels
title('Total Deposition per Year (kg)');
xlabel('Longitude');
ylabel('Latitude');

% Highlight impacted cells
contour(X, Y, impacted_cells', [1 1], 'LineColor', 'r', 'LineWidth', 2);
% Display total deposition
text(lon_max_impacted - offset_x_deg, lat_max_impacted - 1.5*offset_y_deg, ...
    sprintf('Total Deposition:\n%.0f kg', total_deposition_impacted_kg), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'FontWeight', 'bold', 'BackgroundColor', 'w');

% Subplot 6: Population Data (Interpolated)
subplot(1, 2, 2);
hMap3 = imshow(A, 'XData', [lon_min lon_max], 'YData', [lat_min lat_max]);  
set(hMap3, 'AlphaData', 0.5);

% Ensure that the aspect ratio is correct
axis on;        % Turn on axis
axis equal;     % Keep the geographic aspect ratio

hold on;  % Hold on to overlay the contour plot

% Plot the scale bar
plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', 2);
plot([x_start, x_start], [y_start - offset_y_deg * 0.1, y_start + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
plot([x_end, x_end], [y_end - offset_y_deg * 0.1, y_end + offset_y_deg * 0.1], 'k-', 'LineWidth', 2);
text(x_start + scale_bar_length_deg/2, y_start - 0.5*offset_y_deg, ...
    [num2str(scale_length_km) ' km'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'FontWeight', 'bold');

% Create the contour plot for interpolated population data over CFD grid
[c3, cont3] = contourf(X, Y, log10(population_on_CFD_grid), 'LineColor', 'none');
set(gca, 'XLim', [lon_min_impacted lon_max_impacted], 'YLim', [lat_min_impacted lat_max_impacted]);
colorbar;
uistack(cont3, 'bottom');

% Titles and labels
title('Population Count (10³)');
xlabel('Longitude');
ylabel('Latitude');

% Highlight impacted cells on the population map
contour(X, Y, impacted_cells', [1 1], 'LineColor', 'r', 'LineWidth', 2);
% Display total impacted population
text(lon_max_impacted - offset_x_deg, lat_max_impacted - 1.5*offset_y_deg, ...
    sprintf('Total Impacted Population:\n%.0f people', total_population_impacted), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'FontWeight', 'bold', 'BackgroundColor', 'w')