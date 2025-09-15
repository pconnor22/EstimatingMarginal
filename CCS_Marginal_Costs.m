%% This script Preforms the calculation of the average costs for each plant-storage site pair
% Results are in marginal costs for each fossil plant-storage site pair
% (%/metric ton)

E_f = .8./plant_cap(:,2);
plant_cap(:,2) = ones(length(plant_cap),1) .* 0.80;

plant_CO2 = plant_CO2 .* E_f;
plant_CO2_post_capture = plant_CO2_post_capture .* captureRate;

cap_dpt = ones(num_plants,num_sites);

trans_dpt_nom = zeros(num_plants,num_sites);
trans_dpt_PV3 = zeros(num_plants,num_sites);
trans_dpt_PV2 = zeros(num_plants,num_sites);

storage_dpt_nom = zeros(num_plants,num_sites);
storage_dpt_PV3 = zeros(num_plants,num_sites);
storage_dpt_PV2 = zeros(num_plants,num_sites);

storage_CO2 = zeros(num_plants,num_sites);

 [~,c_dpt]= captureModule(plant_CO2, plant_cap(:,2), plant_cap(:,1), plant_fuel,plant_tech,plant_SO2,pload,captureRate,plant_NOX);
 cap_dpt = c_dpt .* cap_dpt;

 clear("c_dpt")

 for i = 1:num_plants
     plant = i
     tstart = tic;

     for j = 1:num_sites

         [~,t_dpt]=transportModule(plant_CO2_post_capture(i),distInMi(i,j)*tortuosity,plant_cap(i,2),escalation_rate,discount_rates);
         trans_dpt_nom(i,j) = t_dpt(1,1);
         trans_dpt_PV3(i,j) = t_dpt(2,1);
         trans_dpt_PV2(i,j) = t_dpt(3,1);

         [~,~,storage_CO2(i,j),s_dpt]=storageModule(PlumeUnArmult,plant_CO2_post_capture(i)*captureRate,j,escalation_rate,discount_rates);
         storage_dpt_nom(i,j) = s_dpt(1,1);
         storage_dpt_PV3(i,j) = s_dpt(2,1);
         storage_dpt_PV2(i,j) = s_dpt(3,1);

     end

     elapsed = toc(tstart)

 end

 total_dpt_nom = cap_dpt + trans_dpt_nom + storage_dpt_nom;
 total_dpt_PV3 = cap_dpt + trans_dpt_PV3 + storage_dpt_PV3;
 total_dpt_PV2 = cap_dpt + trans_dpt_PV2 + storage_dpt_PV2;


 transport_dpt = {trans_dpt_nom;trans_dpt_PV3;trans_dpt_PV2};
 storage_dpt = {storage_dpt_nom;storage_dpt_PV3;storage_dpt_PV2};
 
 total_dpt = {total_dpt_nom;total_dpt_PV3;total_dpt_PV2};
    
 clear("t_dpt","s_dpt","trans_dpt_nom","trans_dpt_PV3","trans_dpt_PV2","storage_dpt_nom","storage_dpt_PV3","storage_dpt_PV2")

 writematrix(total_dpt_nom,"CCS_Outputs/CCS_Tables/2022_Final/total_dpt_nom_1b.csv")
 writematrix(total_dpt_PV3,"CCS_Outputs/CCS_Tables/2022_Final/total_dpt_PV3_1b.csv")
 writematrix(total_dpt_PV2,"CCS_Outputs/CCS_Tables/2022_Final/total_dpt_PV2_1b.csv")

 writematrix(cap_dpt,"CCS_Outputs/CCS_Tables/2022_Final/cap_dpt_1b.csv")

 writematrix(transport_dpt{1,1},"CCS_Outputs/CCS_Tables/2022_Final/trans_dpt_nom_1b.csv")
 writematrix(transport_dpt{2,1},"CCS_Outputs/CCS_Tables/2022_Final/trans_dpt_PV3_1b.csv")
 writematrix(transport_dpt{3,1},"CCS_Outputs/CCS_Tables/2022_Final/trans_dpt_PV2_1b.csv")

 writematrix(storage_dpt{1,1},"CCS_Outputs/CCS_Tables/2022_Final/store_dpt_nom_1b.csv")
 writematrix(storage_dpt{2,1},"CCS_Outputs/CCS_Tables/2022_Final/store_dpt_PV3_1b.csv")
 writematrix(storage_dpt{3,1},"CCS_Outputs/CCS_Tables/2022_Final/store_dpt_PV2_1b.csv")


