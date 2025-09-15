%% Load Plant and Storage Site Data
%This script loads in the plant and storage site location data, along with calculating the distances between 
%the plants  and sites

load("Plant_Site_Data.mat")
pl_lat = plant_latlong(:,1);
pl_lon = plant_latlong(:,2);

% Assuming plant_latlong and latlongSite are matrices of coordinates
earthRadiusInMeters = 6371000; % Earth's radius in meters

num_plants = size(plant_latlong, 1);
num_sites = size(latlongSite, 1);

% Replicate matrices to perform element-wise distance calculations
plant_latlong_rep = repmat(plant_latlong, 1, 1, num_sites);
latlongSite_rep = permute(repmat(latlongSite, 1, 1, num_plants), [3, 2, 1]);

% Compute distances in kilometers
distInKilometers = distance(plant_latlong_rep(:,1,:), plant_latlong_rep(:,2,:), latlongSite_rep(:,1,:), latlongSite_rep(:,2,:), earthRadiusInMeters)/1000;

% Reshape distances  
distInKilometers = reshape(distInKilometers, num_plants, num_sites);

% convert to miles
distInMi = distInKilometers * 0.621371;

plant_CO2_post_capture = plant_CO2.* captureRate;



clear("plant_latlong_rep","latlongSite_rep","plant_latlong","distInKilometers")