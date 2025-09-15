%% Script to take cleaned up data file from R and create a MATLAB workspace
clc
clear

plant_data = readtable('fossil_plants_2022.xlsx');

%Extracting CO2 emissions 
plant_CO2 = table2array(plant_data(:,"PLCO2AN"));

%Extracting plant names
plant_names = string(table2array(plant_data(:,"PNAME")));

%Extracting plant technology - Supercritical, Ultrasupercritical
plant_tech = string(table2array(plant_data(:,47:48)));

%Extracting Plant fuel
plant_fuel = string(table2array(plant_data(:,"PLFUELCT")));

%Extracting plant nameplate capacity & capacity factor
plant_cap = table2array(plant_data(:,29:-1:28));

%Extracting SO2 control
plant_SO2 = string(table2array(plant_data(:,"SO2_Control")));

%Extracting NOX control
plant_NOX = string(table2array(plant_data(:,"NOX_Control")));

%Extracting Age & Retirement
plant_age = table2array(plant_data(:,51:52));

%Extracting Lat & Long
plant_latlong = table2array(plant_data(:,20:21));

%Extracting NOX & SOX
plant_NOX_SOX(:,1) = table2array(plant_data(:,"PLNOXAN"));
plant_NOX_SOX(:,2) = table2array(plant_data(:,"PLSO2AN"));

plant_state = string(table2array(plant_data(:,"PSTATABB")));


%Extraccting ORISPL
plant_ORISPL = double(string(table2array(plant_data(:,"ORISPL"))));

siteDB = readtable('GeolDB.csv'); %Loading in the storage site database

site_state_loc = string(table2array(siteDB(:,4)));

site_basin = string(table2array(siteDB(:,6)));

site_names =string(readcell('sites.txt'));

site_capacity = readtable('capacity.csv');
site_capacity= table2array(site_capacity(table2array(site_capacity(:,1)) == "Reg_dip",2))*1000;
capacity_database = [site_names,site_capacity];

latlongSite = [table2array(siteDB(:,16)),table2array(siteDB(:,15))]; %location of the storage sites

clear plant_data siteDB

save("/Users/patrickconnor/Documents/MATLAB/Graduate School/Research/CCS Model V3/Plant_Site_Data")
