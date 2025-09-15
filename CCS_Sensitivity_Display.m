close all

%% Capture Sensitivity Plots 

c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV3.csv");
st = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV3.csv");
tot = c90 + tr + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_90=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));


c95 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/cap_dpt_95.txt");
tr = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_95.txt");
st = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_95.txt");
tot = c95 + tr + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.95./MMmt);
abate_95=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));

c99 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/cap_dpt_99.txt");
tr = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_99.txt");
st = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_99.txt");
tot = c99 + tr + st;
tot = tot(R_old,:);
[dpt_array,loc] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.99./MMmt);
abate_99=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));


% macTable_scc_abate_99=sortrows(table(dpt_array,plant_CO2_post_capture_array,plant_names));
macTable_scc_abate_99=sortrows(table(dpt_array,plant_CO2_post_capture_array,plant_names,plant_cap(:,1),plant_cap(:,2),plant_ORISPL,plant_fuel,plant_NOX));
writetable(macTable_scc_abate_99,"CCS_Outputs/CCS_MAC/MAC_scc_99.xlsx")

%---------------------------------------------------------------------------

c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV2.csv");
st = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV2.csv");
tot = c90 + tr + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_90_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));


c95 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/cap_dpt_95.txt");
tr = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_95.txt");
st = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_95.txt");
tot = c95 + tr + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.95./MMmt);
abate_95_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));

c99 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/cap_dpt_99.txt");
tr = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_99.txt");
st = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_99.txt");
tot = c99 + tr + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.99./MMmt);
abate_99_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));

macTable_scc_abate_99_2=sortrows(table(dpt_array,plant_CO2_post_capture_array,plant_names,plant_cap(:,1),plant_cap(:,2),plant_ORISPL,plant_fuel,plant_NOX));
writetable(macTable_scc_abate_99_2,"CCS_Outputs/CCS_MAC/MAC_scc_99_2.xlsx")



cap = categorical({'90%','95%','99%'});
cap = reordercats(cap,{'90%','95%','99%'});
y = [abate_90 abate_95 abate_99];

c=barh(cap,y,'black');
c.EdgeColor = "flat";
ylabel("Sensitivity Parameters")
xlabel("Median $/tonne of CO2 across Plants")
c.CData(1,:) = [1 0 0];
c.LineWidth = 3;
saveas(gcf,"CCS_Outputs/CCS_Sensitivity_Figures/Cap_Sensitivity_22.png")
hold on

%% Transport 
c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr10 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_10.txt");
st = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV3.csv");
tot = c90 + tr10 + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_10=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));


c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr15 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV3.csv");
st = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV3.csv");
tot = c90 + tr15 + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_15=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));


c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr20 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV3_20.txt");
st = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV3.csv");
tot = c90 + tr20 + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_20=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));

%----------------------------------------------------------------------------------------

c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr10 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_10.txt");
st = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV2.csv");
tot = c90 + tr10 + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_10_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));


c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr15 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV2.csv");
st = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV2.csv");
tot = c90 + tr15 + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_15_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));


c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr20 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/trans_dpt_PV2_20.txt");
st = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV2.csv");
tot = c90 + tr20 + st;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_20_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));

PV2 = categorical({'10%','15%','20%'});
PV2_t = reordercats(PV2,{'10%','15%','20%'});

PV3 = categorical({'10%','15%','20%'});
PV3 = reordercats(PV3,{'10%','15%',' 20%'});

% bar(PV2,[med_t10_PV2,med_t15_PV2,med_t20_PV2])
hold on 
c=barh(PV3,[abate_10,abate_15,abate_20],'blue');
c.EdgeColor = "flat";
c.CData(2,:) = [1 0 0];
c.LineWidth = 3;
% xlabel("Tortuosity Factor")
% ylabel("Median $/tonne of CO2 across Plants & Storage Sites")
saveas(gcf,"CCS_Outputs/CCS_Sensitivity_Figures/Transport_Sensitivity_22.png")

hold off

%% Storage
c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr15 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV3.csv");
st_175 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV3.csv");
tot = c90 + tr15 + st_175;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_175=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));


c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr15 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV3.csv");
st_1 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_1.txt");
tot = c90 + tr15 + st_1;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_1=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));


c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr15 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV3.csv");
st_225 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV3_225.txt");
tot = c90 + tr15 + st_225;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_225=sum(plant_CO2_post_capture_array(dpt_array <= scc_low));

%----------------------------------------------------------------------------------------

c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr15 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV2.csv");
st_175 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV2.csv");
tot = c90 + tr15 + st_175;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_175_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));


c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr15 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV2.csv");
st_1 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_1.txt");
tot = c90 + tr15 + st_1;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_1_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));


c90 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv");
tr15 = readmatrix("CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV2.csv");
st_225 = readmatrix("CCS_Outputs/CCS_Tables/Sensitivity/2022_SCR/storage_dpt_PV2_225.txt");
tot = c90 + tr15 + st_225;
tot = tot(R_old,:);
[dpt_array,~] = min(tot,[],2);
plant_CO2_post_capture_array=ceil(plant_CO2.*.90./MMmt);
abate_225_2=sum(plant_CO2_post_capture_array(dpt_array <= scc_high));



PV2 = categorical({'1.00','1.75','2.25'});
PV2 = reordercats(PV2,{'1.00','1.75','2.25'});

PV3 = categorical({'1.00','1.75','2.25'});
PV3 = reordercats(PV3,{'1.00','1.75','2.25'});


% bar(PV2,[med_s1_PV2,med_s175_PV2,med_s225_PV2],'yellow')
hold on 
c=barh(PV3,[abate_1,abate_175,abate_225]);
c.EdgeColor = "flat";
c.CData(2,:) = [1 0 0];
c.LineWidth = 3;
xlim([0 275])
% xlabel("Plume Uncertainty Multiplier")
xlabel("Abatement (MMmt)")
legend("Capture Rate","Tortuosity Factor","Plume Uncertainty Multiplier",'Location','northwest');
saveas(gcf,"CCS_Outputs/CCS_Sensitivity_Figures/Storage_Sensitivity_22.png")
hold off


y = [abate_90_2 abate_95_2 abate_99_2];
c=barh(cap,y,'black');
c.EdgeColor = "flat";
c.CData(1,:) = [1 0 0];
c.LineWidth = 3;
xlabel("Sensitivity Parameters")
ylabel("Median $/tonne of CO2 across Plants")
saveas(gcf,"CCS_Outputs/CCS_Sensitivity_Figures/Cap_Sensitivity.png")
hold on
c=barh(PV2_t,[abate_10_2,abate_15_2,abate_20_2],'blue');
c.EdgeColor = "flat";
c.CData(2,:) = [1 0 0];
c.LineWidth = 3;
c=barh(PV2,[abate_1_2,abate_175_2,abate_225_2],'yellow');
c.EdgeColor = "flat";
c.CData(2,:) = [1 0 0];
c.LineWidth = 3;
xlabel("Abatement (MMmt)")
legend("Capture Rate","Tortuosity Factor","Plume Uncertainty Multiplier",'Location','northwest')
ylabel("Sensitivity Parameters")
hold off
saveas(gcf,"CCS_Outputs/CCS_Sensitivity_Figures/Alt Sensitivity 2 per.png")
