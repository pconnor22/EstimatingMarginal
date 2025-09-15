%% This script runs the Integrated Carbon Capture, Transport, and Storage model
% Results are in marginal costs for each fossil plant-storage site pair
% ($/metric ton)

run load_general_assumptions.m
run load_plant_site_data.m

%Creates Cost tables - Takes about 4 hours to run
% run CCS_Marginal_Costs

%Runs CCS_Marginal_Costs w/ Different Sensitivities
% run CCS_Sensitivity.m 

run CCS_Retirement_Old
run CCS_Sensitivity_Display 

% % %MACs
run CCS_MAC_Curve_General
run CCS_MAC_Curve_SCCs
run CCS_MAC_Curve_SCCs_45Q
run CCS_Abate_Basins
% 
% % %Cleanup
run CCS_Cleanup


%Air Quality
% run A_CCS
% addpath("CCS_Air_Quality Work")
% run CCS_Control

%% end of script
