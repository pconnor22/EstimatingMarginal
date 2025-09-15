%% Removing Retiring & Old plants

R_old=plant_age(:,1) <= 50 & ~(plant_age(:,2) < 2030 & ~isnan(plant_age(:,2)));

% R_old=plant_age(:,1) >= 50 | ~(plant_age(:,2) < 2030 & ~isnan(plant_age(:,2)));

% R_old= ~(plant_age(:,2) < 2030 & ~isnan(plant_age(:,2)));

plant_age=plant_age(R_old,:);

plant_cap=plant_cap(R_old,:);

num_plants=length(R_old);

plant_CO2=plant_CO2(R_old);

plant_CO2_post_capture=plant_CO2_post_capture(R_old);

plant_fuel=plant_fuel(R_old);

% plant_latlong=plant_latlong(R_old);

plant_names=plant_names(R_old);

plant_state = plant_state(R_old);

plant_SO2=plant_SO2(R_old);

plant_NOX =plant_NOX(R_old);

plant_tech=plant_tech(R_old,:);

plant_ORISPL=plant_ORISPL(R_old);

plant_NOX_SOX=plant_NOX_SOX(R_old,:);

pl_lat = pl_lat(R_old);

pl_lon = pl_lon(R_old);

