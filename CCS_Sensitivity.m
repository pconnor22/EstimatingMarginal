%% This script preforms the calculation of the average costs for each plant-storage site pair given different sensitive inputs
% Results are in marginal costs for each fossil plant-storage site pair
% (%/metric ton)

%% Capture Rate - 95
disp("CAP - 95")
captureRate = 0.95;

cap_dpt_95 = ones(num_plants,num_sites);

trans_dpt_nom_95 = zeros(num_plants,num_sites);
trans_dpt_PV3_95 = zeros(num_plants,num_sites);
trans_dpt_PV2_95 = zeros(num_plants,num_sites);

storage_dpt_nom_95 = zeros(num_plants,num_sites);
storage_dpt_PV3_95 = zeros(num_plants,num_sites);
storage_dpt_PV2_95 = zeros(num_plants,num_sites);

storage_CO2_95 = zeros(num_plants,num_sites);

[~,c_dpt]= captureModule(plant_CO2, plant_cap(:,2), plant_cap(:,1), plant_fuel,plant_tech,plant_SO2,pload,captureRate,plant_NOX);
cap_dpt_95 = c_dpt .* cap_dpt_95;

for i = 1:num_plants
    plant = i
    tstart = tic;

    for j = 1:num_sites

        [~,t_dpt]=transportModule(plant_CO2(i).*captureRate,distInMi(i,j)*tortuosity,plant_cap(i,2),escalation_rate,discount_rates);
        trans_dpt_nom_95(i,j) = t_dpt(1,1);
        trans_dpt_PV3_95(i,j) = t_dpt(2,1);
        trans_dpt_PV2_95(i,j) = t_dpt(3,1);

        [~,~,storage_CO2_95(i,j),s_dpt]=storageModule(PlumeUnArmult,plant_CO2(i).*captureRate,j,escalation_rate,discount_rates);
        storage_dpt_nom_95(i,j) = s_dpt(1,1);
        storage_dpt_PV3_95(i,j) = s_dpt(2,1);
        storage_dpt_PV2_95(i,j) = s_dpt(3,1);

    end

    elapsed = toc(tstart)

end

total_dpt_nom_95 = cap_dpt_95 + trans_dpt_nom_95 + storage_dpt_nom_95;
total_dpt_PV3_95 = cap_dpt_95 + trans_dpt_PV3_95 + storage_dpt_PV3_95;
total_dpt_PV2_95 = cap_dpt_95 + trans_dpt_PV2_95 + storage_dpt_PV2_95;


writematrix(total_dpt_nom_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/total_dpt_nom_95")
writematrix(total_dpt_PV3_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/total_dpt_PV3_95")
writematrix(total_dpt_PV2_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/total_dpt_PV2_95")

writematrix(cap_dpt_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/cap_dpt_95")
writematrix(storage_CO2_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_CO2_95")


writematrix(storage_dpt_nom_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_nom_95")
writematrix(storage_dpt_PV3_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_95")
writematrix(storage_dpt_PV2_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_95")

writematrix(trans_dpt_nom_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_nom_95")
writematrix(trans_dpt_PV3_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_95")
writematrix(trans_dpt_PV2_95,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_95")

%%


%% This script preforms the calculation of the average costs for each plant-storage site pair given different sensitive inputs
% Results are in marginal costs for each fossil plant-storage site pair
% (%/metric ton)

%% Capture Rate - 99

disp("CAP - 99")
captureRate = 0.99;

cap_dpt_99 = ones(num_plants,num_sites);

trans_dpt_nom_99 = zeros(num_plants,num_sites);
trans_dpt_PV3_99 = zeros(num_plants,num_sites);
trans_dpt_PV2_99 = zeros(num_plants,num_sites);

storage_dpt_nom_99 = zeros(num_plants,num_sites);
storage_dpt_PV3_99 = zeros(num_plants,num_sites);
storage_dpt_PV2_99 = zeros(num_plants,num_sites);

storage_CO2_99 = zeros(num_plants,num_sites);

[~,c_dpt]= captureModule(plant_CO2, plant_cap(:,2), plant_cap(:,1), plant_fuel,plant_tech,plant_SO2,pload,captureRate,plant_NOX);
cap_dpt_99 = c_dpt .* cap_dpt_99;


for i = 1:num_plants
    plant = i
    tstart = tic;

    for j = 1:num_sites

        [~,t_dpt]=transportModule(plant_CO2(i).*captureRate,distInMi(i,j)*tortuosity,plant_cap(i,2),escalation_rate,discount_rates);
        trans_dpt_nom_99(i,j) = t_dpt(1,1);
        trans_dpt_PV3_99(i,j) = t_dpt(2,1);
        trans_dpt_PV2_99(i,j) = t_dpt(3,1);

        [~,~,storage_CO2_99(i,j),s_dpt]=storageModule(PlumeUnArmult,plant_CO2(i).*captureRate,j,escalation_rate,discount_rates);
        storage_dpt_nom_99(i,j) = s_dpt(1,1);
        storage_dpt_PV3_99(i,j) = s_dpt(2,1);
        storage_dpt_PV2_99(i,j) = s_dpt(3,1);

    end

    elapsed = toc(tstart)

end

total_dpt_nom_99 = cap_dpt_99 + trans_dpt_nom_99 + storage_dpt_nom_99;
total_dpt_PV3_99 = cap_dpt_99 + trans_dpt_PV3_99 + storage_dpt_PV3_99;
total_dpt_PV2_99 = cap_dpt_99 + trans_dpt_PV2_99 + storage_dpt_PV2_99;


transport_dpt = {trans_dpt_nom_99;trans_dpt_PV3_99;trans_dpt_PV2_99};
storage_dpt = {storage_dpt_nom_99;storage_dpt_PV3_99;storage_dpt_PV2_99};

total_dpt = {total_dpt_nom_99;total_dpt_PV3_99;total_dpt_PV2_99};


writematrix(total_dpt_nom_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/total_dpt_nom_99")
writematrix(total_dpt_PV3_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/total_dpt_PV3_99")
writematrix(total_dpt_PV2_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/total_dpt_PV2_99")

writematrix(cap_dpt_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/cap_dpt_99")
writematrix(storage_CO2_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_CO2_99")

writematrix(storage_dpt_nom_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_nom_99")
writematrix(storage_dpt_PV3_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_99")
writematrix(storage_dpt_PV2_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_99")

writematrix(trans_dpt_nom_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_nom_99")
writematrix(trans_dpt_PV3_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_99")
writematrix(trans_dpt_PV2_99,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_99")

%% Transport - 10
captureRate = 0.90;

disp("tort - 10")
tortuosity = 1.10;

trans_dpt_nom_10 = zeros(num_plants,num_sites);
trans_dpt_PV3_10 = zeros(num_plants,num_sites);
trans_dpt_PV2_10 = zeros(num_plants,num_sites);


for i = 1:num_plants
    plant = i
    tstart = tic;

    for j = 1:num_sites

        [~,t_dpt]=transportModule(plant_CO2_post_capture(i),distInMi(i,j)*tortuosity,plant_cap(i,2),escalation_rate,discount_rates);
        trans_dpt_nom_10(i,j) = t_dpt(1,1);
        trans_dpt_PV3_10(i,j) = t_dpt(2,1);
        trans_dpt_PV2_10(i,j) = t_dpt(3,1);

    end

    elapsed = toc(tstart)

end

writematrix(trans_dpt_nom_10,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_nom_10")
writematrix(trans_dpt_PV3_10,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_10")
writematrix(trans_dpt_PV2_10,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_10")

%% Transport - 20
disp("tort - 20")
tortuosity = 1.20;

trans_dpt_nom_20 = zeros(num_plants,num_sites);
trans_dpt_PV3_20 = zeros(num_plants,num_sites);
trans_dpt_PV2_20 = zeros(num_plants,num_sites);


for i = 1:num_plants
    plant = i
    tstart = tic;

    for j = 1:num_sites

        [~,t_dpt]=transportModule(plant_CO2_post_capture(i),distInMi(i,j)*tortuosity,plant_cap(i,2),escalation_rate,discount_rates);
        trans_dpt_nom_20(i,j) = t_dpt(1,1);
        trans_dpt_PV3_20(i,j) = t_dpt(2,1);
        trans_dpt_PV2_20(i,j) = t_dpt(3,1);

    end

    elapsed = toc(tstart)

end

writematrix(trans_dpt_nom_20,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_nom_20")
writematrix(trans_dpt_PV3_20,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_20")
writematrix(trans_dpt_PV2_20,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_20")


%% Transport - 30
disp("tort - 30")
tortuosity = 1.30;

trans_dpt_nom_30 = zeros(num_plants,num_sites);
trans_dpt_PV3_30 = zeros(num_plants,num_sites);
trans_dpt_PV2_30 = zeros(num_plants,num_sites);


for i = 1:num_plants
    plant = i
    tstart = tic;

    for j = 1:num_sites

        [~,t_dpt]=transportModule(plant_CO2_post_capture(i),distInMi(i,j)*tortuosity,plant_cap(i,2),escalation_rate,discount_rates);
        trans_dpt_nom_30(i,j) = t_dpt(1,1);
        trans_dpt_PV3_30(i,j) = t_dpt(2,1);
        trans_dpt_PV2_30(i,j) = t_dpt(3,1);

    end

    elapsed = toc(tstart)

end

writematrix(trans_dpt_nom_30,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_nom_30")
writematrix(trans_dpt_PV3_30,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_30")
writematrix(trans_dpt_PV2_30,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_30")

%% store - 1.00
tortuosity=1.15;

PlumeUnArmult = 1;
disp("storage - 1")

storage_dpt_nom_1 = zeros(num_plants,num_sites);
storage_dpt_PV3_1 = zeros(num_plants,num_sites);
storage_dpt_PV2_1 = zeros(num_plants,num_sites);

storage_CO2_1 = zeros(num_plants,num_sites);


 for i = 1:num_plants
     plant = i
     tstart = tic;

     for j = 1:num_sites

         [~,~,storage_CO2_1(i,j),s_dpt]=storageModule(PlumeUnArmult,plant_CO2_post_capture(i)*captureRate,j,escalation_rate,discount_rates);
         storage_dpt_nom_1(i,j) = s_dpt(1,1);
         storage_dpt_PV3_1(i,j) = s_dpt(2,1);
         storage_dpt_PV2_1(i,j) = s_dpt(3,1);

     end

     elapsed = toc(tstart)

 end

 writematrix(storage_CO2_1,"CCS_Outputs/CCS_Tables/Sensitivity/storage_CO2_1")

 writematrix(storage_dpt_nom_1,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_1")
 writematrix(storage_dpt_PV3_1,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_1")
 writematrix(storage_dpt_PV2_1,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_1")

% store - 2.25
tortuosity=1.15;

PlumeUnArmult = 2.25;
disp("storage - 225")

storage_dpt_nom_225 = zeros(num_plants,num_sites);
storage_dpt_PV3_225 = zeros(num_plants,num_sites);
storage_dpt_PV2_225 = zeros(num_plants,num_sites);

storage_CO2_225 = zeros(num_plants,num_sites);


 for i = 1:num_plants
     plant = i
     tstart = tic;

     for j = 1:num_sites

         [~,~,storage_CO2_225(i,j),s_dpt]=storageModule(PlumeUnArmult,plant_CO2_post_capture(i)*captureRate,j,escalation_rate,discount_rates);
         storage_dpt_nom_225(i,j) = s_dpt(1,1);
         storage_dpt_PV3_225(i,j) = s_dpt(2,1);
         storage_dpt_PV2_225(i,j) = s_dpt(3,1);

     end

     elapsed = toc(tstart)

 end

 writematrix(storage_CO2_225,"CCS_Outputs/CCS_Tables/Sensitivity/storage_CO2_225")

 writematrix(storage_dpt_nom_225,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_nom_225")
 writematrix(storage_dpt_PV3_225,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_225")
 writematrix(storage_dpt_PV2_225,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_225")


 %% store - 1.75
PlumeUnArmult = 1.75;
disp("storage - 175")

storage_dpt_nom_175 = zeros(num_plants,num_sites);
storage_dpt_PV3_175 = zeros(num_plants,num_sites);
storage_dpt_PV2_175 = zeros(num_plants,num_sites);

storage_CO2_225 = zeros(num_plants,num_sites);


 for i = 1:num_plants
     plant = i
     tstart = tic;

     for j = 1:num_sites

         [~,~,storage_CO2_175(i,j),s_dpt]=storageModule(PlumeUnArmult,plant_CO2_post_capture(i)*captureRate,j,escalation_rate,discount_rates);
         storage_dpt_nom_175(i,j) = s_dpt(1,1);
         storage_dpt_PV3_175(i,j) = s_dpt(2,1);
         storage_dpt_PV2_175(i,j) = s_dpt(3,1);

     end

     elapsed = toc(tstart)

 end

 writematrix(storage_CO2_175,"CCS_Outputs/CCS_Tables/Sensitivity/storage_CO2_175")

 writematrix(storage_dpt_nom_175,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_nom_175")
 writematrix(storage_dpt_PV3_175,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_175")
 writematrix(storage_dpt_PV2_175,"CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_175")






