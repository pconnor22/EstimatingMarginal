%% Script to create Capture Workspace

clear
clc

learning = "Y";
operation_yrs = 30;

lr_OM_NGCC = 0.039;
lr_OM_PC = 0.057;
% lr_factor = floor(log2(operation_yrs));


%% Subcritical Coal
%Newly built greenfield plant data
sub_coal_wo_capture_toc = 1614635;
sub_coal_wo_capture_fixed_OM = 44356511;
sub_coal_wo_capture_variable_OM = 38157907;
sub_coal_wo_capture_fuel_use = 5917;

%Retrofit Plant Data
%Additional capture capital costs
%Costs were chosen where there were extreme differences and/or unique to
%the capture case
FGD_sub = 134771+3137;
FGD_sub_VOC = 3906999 + 2425548;
SCR_sub = 66405;
SCR_sub_VOC = 774077 + 12901;


CO2_removal_system = 544602;
CO2_compression_drying = 91673;
CO2_compressor_aftercooler = 1085;
sub_coal_w_capture_toc = CO2_removal_system + CO2_compression_drying + CO2_compressor_aftercooler;

CO2_system_maintenance_labor = 8384879;
sub_coal_w_capture_fixed_OM = sub_coal_wo_capture_fixed_OM + CO2_system_maintenance_labor;

CO2_maintenance_material = 12577319;
CO2_chemicals = 8857728;
Triethylene_glycol_use = 1321278;
Triethylene_glycol_disposal = 68007;
Thermal_reclaimer_unit_waste = 95098;
sub_coal_w_capture_variable_OM = sub_coal_wo_capture_variable_OM + CO2_maintenance_material + CO2_chemicals + Triethylene_glycol_use +...
    Triethylene_glycol_disposal + Thermal_reclaimer_unit_waste;

sub_coal_w_capture_fuel_use = sub_coal_wo_capture_fuel_use;

%% Supercritical coal
sup_coal_wo_capture_toc = 1682505;
sup_coal_wo_capture_fixed_OM = 45959058;
sup_coal_wo_capture_variable_OM = 37472784;
sup_coal_wo_capture_fuel_use = 5677;

FGD_sup = 130040+3021;
FGD_sup_VOC = 4760201 + 2955235;
SCR_sup = 64193;
SCR_sup_VOC = 737100 + 12285;

CO2_removal_system = 528109;
CO2_compression_drying = 88962;
CO2_compressor_aftercooler = 1045;
sup_coal_w_capture_toc = CO2_removal_system + CO2_compression_drying + CO2_compressor_aftercooler;

CO2_system_maintenance_labor = 8751176;
sup_coal_w_capture_fixed_OM = sup_coal_wo_capture_fixed_OM + CO2_system_maintenance_labor;

CO2_maintenance_material = 13126763;
CO2_chemicals = 8434609;
Triethylene_glycol_use = 1258163;
Triethylene_glycol_disposal = 64758;
Thermal_reclaimer_unit_waste = 90555;
sup_coal_w_capture_variable_OM = sup_coal_wo_capture_variable_OM + CO2_maintenance_material + CO2_chemicals + Triethylene_glycol_use +...
    Triethylene_glycol_disposal + Thermal_reclaimer_unit_waste;

sup_coal_w_capture_fuel_use = sup_coal_wo_capture_fuel_use;

%% Gas Plants
NG_wo_capture_toc = 691670;
NG_wo_capture_fixed_OM = 19466549;
NG_wo_capture_variable_OM = 9225789;
NG_wo_capture_fuel_use = 110955;

SCR_gas = 4057;
SCR_gas_VOC = 144041 + 2401;

CO2_removal_system = 415297;
CO2_compression_drying = 60817;
CO2_compressor_aftercooler = 513;
Gas_cleanup_foundation = 1145;
NG_w_capture_toc = CO2_removal_system + CO2_compression_drying + CO2_compressor_aftercooler + Gas_cleanup_foundation;

CO2_system_maintenance_labor = 4309151;
NG_w_capture_fixed_OM = NG_wo_capture_fixed_OM + CO2_system_maintenance_labor;

CO2_maintenance_material = 6463727;
CO2_chemicals = 1969571;
Triethylene_glycol_use = 931690;
Triethylene_glycol_disposal = 47955;
Thermal_reclaimer_unit_waste = 20643;
NG_w_capture_variable_OM = NG_wo_capture_variable_OM + CO2_maintenance_material + CO2_chemicals + Triethylene_glycol_use +...
    Triethylene_glycol_disposal + Thermal_reclaimer_unit_waste;

NG_w_capture_fuel_use = NG_wo_capture_fuel_use;


if learning =="Y"
    avg_fixed_OM_coal = 0;
    avg_var_OM_coal = 0;

    avg_fixed_OM_gas = 0;
    avg_var_OM_gas = 0;

    for i=1:operation_yrs
        avg_fixed_OM_coal = avg_fixed_OM_coal+((1-lr_OM_PC)^log2(i)) * (sub_coal_w_capture_fixed_OM - sub_coal_wo_capture_fixed_OM);
        avg_var_OM_coal = avg_var_OM_coal + ((1-lr_OM_PC)^log2(i)) * (sub_coal_w_capture_variable_OM - sub_coal_wo_capture_variable_OM);
        
        avg_fixed_OM_gas = avg_fixed_OM_gas+((1-lr_OM_NGCC)^log2(i)) * (NG_w_capture_fixed_OM - NG_wo_capture_fixed_OM);
        avg_var_OM_gas = avg_var_OM_gas + ((1-lr_OM_NGCC)^log2(i)) * (NG_w_capture_variable_OM - NG_wo_capture_variable_OM);
    end

     avg_fixed_OM_coal = avg_fixed_OM_coal/30;
    avg_var_OM_coal = avg_var_OM_coal/30;

    avg_fixed_OM_gas = avg_fixed_OM_gas/30;
    avg_var_OM_gas = avg_var_OM_gas/30;

    sub_coal_w_capture_fixed_OM = sub_coal_wo_capture_fixed_OM + avg_fixed_OM_coal;
    sub_coal_w_capture_variable_OM = sub_coal_wo_capture_variable_OM + avg_var_OM_coal;

    sup_coal_w_capture_fixed_OM = sup_coal_wo_capture_fixed_OM + avg_fixed_OM_coal;
    sup_coal_w_capture_variable_OM = sup_coal_wo_capture_variable_OM + avg_var_OM_coal;

    NG_w_capture_fixed_OM = NG_wo_capture_fixed_OM + avg_fixed_OM_gas;
    NG_w_capture_variable_OM = NG_wo_capture_variable_OM + avg_var_OM_gas;

end




%% Packaging

sub_coal = [sub_coal_wo_capture_toc,sub_coal_w_capture_toc;
    sub_coal_wo_capture_fixed_OM,sub_coal_w_capture_fixed_OM;
    sub_coal_wo_capture_variable_OM,sub_coal_w_capture_variable_OM;
    sub_coal_wo_capture_fuel_use,sub_coal_w_capture_fuel_use];

sup_coal = [sup_coal_wo_capture_toc,sup_coal_w_capture_toc;
    sup_coal_wo_capture_fixed_OM,sup_coal_w_capture_fixed_OM;
    sup_coal_wo_capture_variable_OM,sup_coal_w_capture_variable_OM;
    sup_coal_wo_capture_fuel_use,sup_coal_w_capture_fuel_use];

NG = [NG_wo_capture_toc,NG_w_capture_toc;
    NG_wo_capture_fixed_OM,NG_w_capture_fixed_OM;
    NG_wo_capture_variable_OM,NG_w_capture_variable_OM;
    NG_wo_capture_fuel_use,NG_w_capture_fuel_use];

FGD_units = [FGD_sub,FGD_sub_VOC,FGD_sup,FGD_sup_VOC];
SCR_units = [SCR_sub,SCR_sub_VOC,SCR_sup,SCR_sup_VOC,SCR_gas,SCR_gas_VOC];

variablesToKeep = {'NG', 'sup_coal', 'sub_coal','FGD_units','SCR_units'};

% Clear all variables except the ones specified
allVars = who;  % Get a list of all variables in the workspace
varsToClear = setdiff(allVars, variablesToKeep);  % Find variables to clear
clear(varsToClear{:});  % Clear variables except those specified
clear("varsToClear","allVars")

%Fixed charged rate, Total as spent cost - total overnight cost ratio, OM
%levelization, coal and gas prices
FCR = 0.088564060;
TASC_TOC_35 = 1.289;
TASC_TOC_33 = 1.242;
OMlevel = 1.384;
illinoisCoal = 2.2269;
midwestNaturalGas = 4.42;



save("/Users/patrickconnor/Documents/MATLAB/Graduate School/Research/CCS Model V3/Capture_Data_12")



