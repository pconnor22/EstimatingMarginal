%% This script takes the outputted cost tables and assembles them into Marginal Abatement Cost Curves by minimizing across storage sites and exaimines the effect of 45Q
%Curves include SCC values & Abatement Levels along with 45Q
close all

%% 12 years of 45Q  3%

% dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022/total_dpt_PV3.csv');
dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/total_dpt_PV3_lr.csv');
% dpt_matrix(~R_old,:) = inf;
dpt_matrix = dpt_matrix(R_old,:);

plt=(1:length(plant_names))';


%Minimizing
[dpt_array,loc] = min(dpt_matrix,[],2);

sub = repmat(dpt_45Q(1,1),length(dpt_array),1);
dpt_array_45Q = dpt_array - sub;

%Converting Emission data into a row vector
plant_CO2_post_capture_array=ceil(plant_CO2_post_capture/MMmt);

sites_array = site_names(loc);
basins = site_basin(loc);
st = site_state_loc(loc);


%turning the 4 row vectors into a table and then sorting
macTable_45Q_12 =sortrows(table(dpt_array,dpt_array_45Q,plant_CO2_post_capture_array,plant_names,plant_fuel,sites_array,basins,st,plant_ORISPL,plt,plant_age(:,1),plant_cap(:,2),loc,pl_lat,pl_lon,plant_NOX,plant_cap(:,1),plant_state));

sortedDPT = table2array(macTable_45Q_12(:,1));
sortedDPT45 = table2array(macTable_45Q_12(:,2));
sortedEmiss = table2array(macTable_45Q_12(:,3));

%creating a cell array to store the $/tonne vectors for the total
%tonnes for each source
out = cell(1,numel(sortedDPT));
out2 = cell(1,numel(sortedDPT45));

tstart = tic;
%looping through each source sink pair
for i = 1:length(sortedDPT)
    %excluding $/tonne costs that are infinite or egregious (> 1e8) due to
    %their distortionary effect on the curve
    if sortedDPT(i) == inf || sortedDPT(i) > 250
    else
    %filling the specific cell of the cell array with a vector that is the
    %length of total tonnes from a source and filled  with marginal
    %cost($/tonne) of each tonne stored
    out{i} = repmat(sortedDPT(i),1,sortedEmiss(i));
    end
    if sortedDPT45(i) == inf || sortedDPT45(i) > 250
    else
    %filling the specific cell of the cell array with a vector that is the
    %length of total tonnes from a source and filled  with marginal
    %cost($/tonne) of each tonne stored
    out2{i} = repmat(sortedDPT45(i),1,sortedEmiss(i));
    end
end

%convertingg the cell arrray into one long double row vector
macCurve45_12 = horzcat(out2{:});
%axis_line = zeros(1,length(macCurve));


%plotting the curve
plot(macCurve_scc_low,'LineWidth',5)%original pv3 maccurve
hold on
plot(macCurve45_12,'LineWidth',5)
% plot(axis_line,'LineWidth',5,'Color','black','LineStyle','--')
xlabel('Abatement (MMmt CO_2)')
ylabel('Marginal Abatement Cost ($/tCO2)')
% legend('MAC-45Q-PV3','Axis Line','Location','northwest')
%title('Marginal Abatement Curve')
saveas(gcf,"CCS_Outputs/CCS_MAC/MAC_12_45Q_PV3.png")
writetable(macTable_45Q_12,"CCS_Outputs/CCS_MAC/MAC_12_45Q_PV3.xlsx")

retrofits = table2array(macTable_45Q_12(:,2));
retrofits = retrofits(retrofits < 0);

num_retro = length(retrofits);
rt_fuel = table2array(macTable_45Q_12(:,5));
rt_gas = sum(rt_fuel(retrofits < 0) == "GAS");
rt_coal = sum(rt_fuel(retrofits < 0) == "COAL");
abate = table2array(macTable_45Q_12(:,3));
c_f =  table2array(macTable_45Q_12(:,12));
pl_age =  table2array(macTable_45Q_12(:,11));
rt_cost = table2array(macTable_45Q_12(:,1));



disp("12 Years 45Q w/ Present Value Costs discounted @ 3%:")
disp("----------------------------------------")
disp(["# Retrofits:",num2str(num_retro)])
disp(["# Coal:",num2str(rt_coal)])
disp(["# Gas:",num2str(rt_gas)])
disp(["Amnt Abated:",num2str(sum(abate(retrofits < 0)))])
disp(["% 2022 Emissions:",num2str(sum(abate(retrofits < 0))/Emiss_2022*100)]);
disp(["Avg CF:",num2str(mean(c_f(retrofits < 0)))]);
disp(["Avg tonnage:",num2str(mean(abate(retrofits < 0)))]);
disp(["Avg Age:",num2str(mean(pl_age(retrofits < 0)))]);
disp(["Median total cost:",num2str(median(abate(retrofits < 0).*MMmt.*rt_cost(retrofits<0)))]);


%Capacity
site_cap_45PV3 = site_capacity;
pl_CO2_45PV3 = table2array(macTable_45Q_12(:,3));
sites_45PV3 = table2array(macTable_45Q_12(:,13));

j=1;
while sum(site_cap_45PV3 < 0) == 0

for i=1:num_retro
    site_cap_45PV3(sites_45PV3(i)) = site_cap_45PV3(sites_45PV3(i)) - pl_CO2_45PV3(i);

end
j=j+1;

end

disp("Capacity Information ")
disp(["It would take ",num2str(j)," years to fill a site"])
filled_sites = find(site_cap_45PV3 < 0);
names_filled_sites= site_names(find(site_cap_45PV3 < 0));
disp(table(filled_sites,names_filled_sites));


%% 45Q - Discount 3%

% dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022/total_dpt_PV3.csv');
dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/total_dpt_PV3_lr.csv');

% plt=(1:length(plant_names))';
plt = 1:length(dpt_matrix);

dpt_matrix = dpt_matrix(R_old,:);
% dpt_matrix(~R_old,:) = inf;

plt=plt(R_old)';


%Minimizing
[dpt_array,loc] = min(dpt_matrix,[],2);

sub = repmat(dpt_45Q(1,2),length(dpt_array),1);
dpt_array_45Q = dpt_array - sub;

%Converting Emission data into a row vector
plant_CO2_post_capture_array=ceil(plant_CO2_post_capture/MMmt);

sites_array = site_names(loc);
basins = site_basin(loc);
st = site_state_loc(loc);


%turning the 4 row vectors into a table and then sorting
macTable_45Q_PV3 =sortrows(table(dpt_array,dpt_array_45Q,plant_CO2_post_capture_array,plant_names,plant_fuel,sites_array,basins,st,plant_ORISPL,plt,plant_age(:,1),plant_cap(:,2),loc,pl_lat,pl_lon,plant_NOX,plant_cap(:,1),plant_state));

sortedDPT = table2array(macTable_45Q_PV3(:,1));
sortedDPT45 = table2array(macTable_45Q_PV3(:,2));
sortedEmiss = table2array(macTable_45Q_PV3(:,3));

%creating a cell array to store the $/tonne vectors for the total
%tonnes for each source
out = cell(1,numel(sortedDPT));
out2 = cell(1,numel(sortedDPT45));

tstart = tic;
%looping through each source sink pair
for i = 1:length(sortedDPT)
    %excluding $/tonne costs that are infinite or egregious (> 1e8) due to
    %their distortionary effect on the curve
    if sortedDPT(i) == inf || sortedDPT(i) > 250
    else
    %filling the specific cell of the cell array with a vector that is the
    %length of total tonnes from a source and filled  with marginal
    %cost($/tonne) of each tonne stored
    out{i} = repmat(sortedDPT(i),1,sortedEmiss(i));
    end
    if sortedDPT45(i) == inf || sortedDPT45(i) > 250
    else
    %filling the specific cell of the cell array with a vector that is the
    %length of total tonnes from a source and filled  with marginal
    %cost($/tonne) of each tonne stored
    out2{i} = repmat(sortedDPT45(i),1,sortedEmiss(i));
    end
end

%convertingg the cell arrray into one long double row vector
macCurve45 = horzcat(out2{:});
axis_line = zeros(1,length(macCurve45));


%plotting the curve
plot(macCurve45,'LineWidth',5)
hold on
plot(axis_line,'LineWidth',5,'Color','black','LineStyle','-.')
xlabel('Abatement (MMmt CO_2)','Interpreter','tex')
ylabel('Marginal Abatement Cost ($/mtCO_2)','Interpreter','tex')
legend('MAC-45Q-PV3','Axis Line','Location','northwest')
legend('MAC','MAC - 45Q 12-year lifetime','MAC - 45Q Full lifetime','Zero-cost Line','Location','northwest')
%title('Marginal Abatement Curve')
% ax = gca;
% ax.XAxis.TickLabelFormat = '%,g';
% ax.YAxis.TickLabelFormat = '%,g';
% axis equal
saveas(gcf,"CCS_Outputs/CCS_MAC/MAC_45Q_PV3.png")
writetable(macTable_45Q_PV3,"CCS_Outputs/CCS_MAC/MAC_45Q_PV3.xlsx")

retrofits = table2array(macTable_45Q_PV3(:,2));
retrofits = retrofits(retrofits < 0);

num_retro = length(retrofits);
rt_fuel = table2array(macTable_45Q_PV3(:,5));
rt_gas = sum(rt_fuel(retrofits < 0) == "GAS");
rt_coal = sum(rt_fuel(retrofits < 0) == "COAL");
abate = table2array(macTable_45Q_PV3(:,3));
c_f =  table2array(macTable_45Q_PV3(:,12));
pl_age =  table2array(macTable_45Q_PV3(:,11));
rt_cost = table2array(macTable_45Q_PV3(:,1));




disp("45Q w/ Present Value Costs discounted @ 3%:")
disp("----------------------------------------")
disp(["# Retrofits:",num2str(num_retro)])
disp(["# Coal:",num2str(rt_coal)])
disp(["# Gas:",num2str(rt_gas)])
disp(["Amnt Abated:",num2str(sum(abate(retrofits < 0)))])
disp(["% 2022 Emissions:",num2str(sum(abate(retrofits < 0))/Emiss_2022*100)]);
disp(["Avg CF:",num2str(mean(c_f(retrofits < 0)))]);
disp(["Avg tonnage:",num2str(mean(abate(retrofits < 0)))]);
disp(["Avg Age:",num2str(mean(pl_age(retrofits < 0)))]);
disp(["Median total cost:",num2str(median(abate(retrofits < 0).*MMmt.*rt_cost(retrofits<0)))]);

%Capacity
site_cap_45PV3 = site_capacity;
pl_CO2_45PV3 = table2array(macTable_45Q_PV3(:,3));
sites_45PV3 = table2array(macTable_45Q_PV3(:,13));

j=1;
while sum(site_cap_45PV3 < 0) == 0

for i=1:num_retro
    site_cap_45PV3(sites_45PV3(i)) = site_cap_45PV3(sites_45PV3(i)) - pl_CO2_45PV3(i);

end
j=j+1;

end

disp("Capacity Information ")
disp(["It would take ",num2str(j)," years to fill a site"])
filled_sites = find(site_cap_45PV3 < 0);
names_filled_sites= site_names(find(site_cap_45PV3 < 0));
disp(table(filled_sites,names_filled_sites));
