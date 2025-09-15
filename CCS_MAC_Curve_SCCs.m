%% This script takes the outputted cost tables and assembles them into Marginal Abatement Cost Curves by minimizing across storage sites
%Curves include SCC values & Abatement Levels
close all

%% SCC - Low

dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/total_dpt_PV3_lr.csv');
plt=(1:length(plant_names))';

cap_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv');
trans_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV3.csv');
store_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV3.csv');

cap_matrix = cap_matrix(R_old,:);
trans_matrix = trans_matrix(R_old,:);
store_matrix = store_matrix(R_old,:);

% dpt_matrix = readmatrix('CCS_Outputs/total_dpt_PV2.csv');
dpt_matrix = dpt_matrix(R_old,:);
% dpt_matrix(~R_old,:) = inf;

% cap_matrix = cap_matrix(R_old,:);
% trans_matrix = trans_matrix(R_old,:);
% store_matrix = store_matrix(R_old,:);

%Minimizing
[dpt_array,loc] = min(dpt_matrix,[],2);
%Converting Emission data into a row vector
plant_CO2_post_capture_array=ceil(plant_CO2_post_capture/MMmt);

sites_array = site_names(loc);
basins = site_basin(loc);
st = site_state_loc(loc);

cap_array = cap_matrix((1:length(loc))' + (loc(:,1)-1)*size(cap_matrix,1));
trans_array = trans_matrix((1:length(loc))' + (loc(:,1)-1)*size(trans_matrix,1));
store_array = store_matrix((1:length(loc))' + (loc(:,1)-1)*size(store_matrix,1));

%turning the 4 row vectors into a table and then sorting
macTable_scc_low =sortrows(table(dpt_array,plant_CO2_post_capture_array,plant_names,plant_fuel,sites_array,basins,st,plant_ORISPL,plt,plant_age(:,1),plant_cap(:,2),loc,cap_array,trans_array,store_array,pl_lat,pl_lon,plant_NOX,plant_cap(:,1),plant_state));

sortedDPT = table2array(macTable_scc_low(:,1));
sortedEmiss = table2array(macTable_scc_low(:,2));

%creating a cell array to store the $/tonne vectors for the total
%tonnes for each source
out = cell(1,numel(sortedDPT));

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
end

%convertingg the cell arrray into one long double row vector
macCurve_scc_low = horzcat(out{:});
SCCnow = repmat(scc_low,1,length(macCurve_scc_low));

%plotting the curve
%subplot(2,1,1);

% figure('Position', [100, 100, 700, 700]); % [left, bottom, width, height]
tiledlayout(2,1)
nexttile;
plot(macCurve_scc_low,'LineWidth',5)
hold on
plot(SCCnow,'LineWidth',5,'Color','black','LineStyle','-.')
% plot(SCCnew,'LineWidth',5,'Color','black','LineStyle','-.')
xlabel('Abatement (MMmt CO_{2})',Interpreter='tex')
% ylabel('\textbf{Marginal Abatement Cost (\$/mt${CO_2}$)}','Interpreter','latex')
ylabel('Marginal Abatement Cost ($/mtCO_{2})','Interpreter','tex')
% legend('MAC',['SCC-',num2str(scc_low)],'Location','northwest')
legend('MAC',['MB-',num2str(scc_low)],'Location','northwest')
%title('Marginal Abatement Curve')
hold off
% axis equal;
% ax = gca;
% ax.XAxis.TickLabelFormat = '%,g';
% ax.YAxis.TickLabelFormat = '%,g';
saveas(gcf,"CCS_Outputs/CCS_MAC/MAC_scc_low.png")
writetable(macTable_scc_low,"CCS_Outputs/CCS_MAC/MAC_scc_low.xlsx")

retrofits = table2array(macTable_scc_low(:,1));
retrofits = retrofits(retrofits < scc_low);

num_retro = length(retrofits);
rt_fuel = table2array(macTable_scc_low(:,4));
rt_gas = sum(rt_fuel(retrofits < scc_low) == "GAS");
rt_coal = sum(rt_fuel(retrofits < scc_low) == "COAL");
abate = table2array(macTable_scc_low(:,2));
c_f =  table2array(macTable_scc_low(:,11));
pl_age =  table2array(macTable_scc_low(:,10));
rt_cost = table2array(macTable_scc_low(:,1));
abate_cost = trapz(macCurve_scc_low(macCurve_scc_low < scc_low))*1e6;
capacity = table2array(macTable_scc_low(:,1));
tonnage = table2array(macTable_scc_low(:,2));


disp(["In an SCC Scenario of",num2str(scc_low)," $/tonne of CO2:"])
disp("----------------------------------------")
disp(["# Retrofits:",num2str(num_retro)])
disp(["# Coal:",num2str(rt_coal)])
disp(["# Gas:",num2str(rt_gas)])
disp(["Amnt Abated:",num2str(sum(abate(retrofits < scc_low)))])
disp(["% 2022 Emissions:",num2str(sum(abate(retrofits < scc_low))/Emiss_2022*100)]);
disp(["Avg CF:",num2str(mean(c_f(retrofits < scc_low)))]);
disp(["Avg tonnage:",num2str(mean(tonnage(retrofits < scc_low)))]);
% disp(["Capacity:", num2str(mean())])
disp(["Avg Age:",num2str(mean(pl_age(retrofits < scc_low)))]);
disp(["Median total cost:",num2str(median(abate(retrofits < scc_low).*MMmt.*rt_cost(retrofits<scc_low)))]);
disp(["Total Abatement Cost ($ Billion):", num2str(abate_cost/1e9)]);

%Capacity
site_cap_scc_low = site_capacity;
pl_CO2_scc_low = table2array(macTable_scc_low(:,2));
sites_scc_low = table2array(macTable_scc_low(:,12));

j=1;
while sum(site_cap_scc_low < 0) == 0

for i=1:num_retro
    site_cap_scc_low(sites_scc_low(i)) = site_cap_scc_low(sites_scc_low(i)) - pl_CO2_scc_low(i);

end
j=j+1;

end


disp("Capacity Information ")
disp(["It would take ",num2str(j)," years to fill a site"])
filled_sites = find(site_cap_scc_low < 0);
names_filled_sites= site_names(find(site_cap_scc_low < 0));
disp(table(filled_sites,names_filled_sites));



%% SCC - High

dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/total_dpt_PV2_lr.csv');
dpt_matrix = dpt_matrix(R_old,:);
% dpt_matrix(~R_old,:) = inf;

cap_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/cap_dpt_lr.csv');
trans_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/trans_dpt_PV2.csv');
store_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/store_dpt_PV2.csv');

cap_matrix = cap_matrix(R_old,:);
trans_matrix = trans_matrix(R_old,:);
store_matrix = store_matrix(R_old,:);


%Minimizing
[dpt_array,loc] = min(dpt_matrix,[],2);

%Converting Emission data into a row vector
plant_CO2_post_capture_array=ceil(plant_CO2_post_capture/MMmt);

sites_array = site_names(loc);
basins = site_basin(loc);
st = site_state_loc(loc);

cap_array = cap_matrix((1:length(loc))' + (loc(:,1)-1)*size(cap_matrix,1));
trans_array = trans_matrix((1:length(loc))' + (loc(:,1)-1)*size(trans_matrix,1));
store_array = store_matrix((1:length(loc))' + (loc(:,1)-1)*size(store_matrix,1));


%turning the 4 row vectors into a table and then sorting
macTable_scc_high =sortrows(table(dpt_array,plant_CO2_post_capture_array,plant_names,plant_fuel,sites_array,basins,st,plant_ORISPL,plt,plant_age(:,1),plant_cap(:,2),loc,cap_array,trans_array,store_array,pl_lat,pl_lon,plant_NOX,plant_cap(:,1),plant_state));

sortedDPT = table2array(macTable_scc_high(:,1));
sortedEmiss = table2array(macTable_scc_high(:,2));

%creating a cell array to store the $/tonne vectors for the total
%tonnes for each source
out = cell(1,numel(sortedDPT));

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
end

%convertingg the cell arrray into one long double row vector
macCurve_scc_high = horzcat(out{:});
SCCnew = repmat(scc_high,1,length(macCurve_scc_high));

%plotting the curve
%subplot(2,1,2);
nexttile;
plot(macCurve_scc_high,'LineWidth',5);
hold on
% plot(SCCnow,'LineWidth',5,'Color','black','LineStyle','--')
plot(SCCnew,'LineWidth',5,'Color','black','LineStyle','-.')
xlabel('Abatement (MMmt CO_{2})')
ylabel('Marginal Abatement Cost ($/mtCO_{2})','Interpreter','tex')
% legend('MAC',['SCC-',num2str(scc_high)],'Location','northwest')
legend('MAC',['MB-',num2str(scc_high)],'Location','northwest')
% axis equal;
% %title('Marginal Abatement Curve')
% ax = gca;
% ax.XAxis.TickLabelFormat = '%,g';
% ax.YAxis.TickLabelFormat = '%,g';

saveas(gcf,"CCS_Outputs/CCS_MAC/MAC_scc_high.png")
%print('CCS_Outputs/CCS_MAC/MAC_scc_high.pdf','-dpdf','-bestfit')%,'resize','-bestfit')
writetable(macTable_scc_high,"CCS_Outputs/CCS_MAC/MAC_scc_high.xlsx")

retrofits = table2array(macTable_scc_high(:,1));
retrofits = retrofits(retrofits < scc_high);

num_retro = length(retrofits);
rt_fuel = table2array(macTable_scc_high(:,4));

gas_loc = rt_fuel(retrofits < scc_high) == "GAS";
coal_loc = rt_fuel(retrofits < scc_high) == "COAL";
rt_gas = sum(gas_loc);
rt_coal = sum(coal_loc);
abate = table2array(macTable_scc_high(:,2));
abate_cost = trapz(macCurve_scc_high(macCurve_scc_high < scc_high))*1e6;


hold off


disp(["In an SCC Scenario of",num2str(scc_high)," $/tonne of CO2:"])
disp("----------------------------------------")
disp(["# Retrofits:",num2str(num_retro)])
disp(["# Coal:",num2str(rt_coal)])
disp(["# Gas:",num2str(rt_gas)])
disp(["Amnt Abated:",num2str(sum(abate(retrofits < scc_high)))])
disp(["% 2022 Emissions:",num2str(sum(abate(retrofits < scc_high))/Emiss_2022*100)]);
disp(["Total Abatement Cost ($ Billion):", num2str(abate_cost/1e9)]);


disp("Gas")
disp(["Amnt Abated:",num2str(sum(abate(gas_loc)))])
disp(["Avg ton:", num2str(mean(abate(gas_loc)))])
disp(["Avg CF:",num2str(mean(c_f(gas_loc)))]);
disp(["Avg Age:",num2str(mean(pl_age(gas_loc)))]);
disp(["Median total cost:",num2str(median(abate(gas_loc).*MMmt.*rt_cost(gas_loc)))]);

disp("Coal")
disp(["Amnt Abated:",num2str(sum(abate(coal_loc)))])
disp(["Avg ton:", num2str(mean(abate(coal_loc)))])
disp(["Avg CF:",num2str(mean(c_f(coal_loc)))]);
disp(["Avg Age:",num2str(mean(pl_age(coal_loc)))]);
disp(["Median total cost:",num2str(median(abate(coal_loc).*MMmt.*rt_cost(coal_loc)))]);




%Capacity
site_cap_scc_high = site_capacity;
pl_CO2_scc_high = table2array(macTable_scc_high(:,2));
sites_scc_high = table2array(macTable_scc_high(:,12));

j=1;
while sum(site_cap_scc_high < 0) == 0

for i=1:num_retro
    site_cap_scc_high(sites_scc_high(i)) = site_cap_scc_high(sites_scc_high(i)) - pl_CO2_scc_high(i);

end
j=j+1;

end

disp("Capacity Information ")
disp(["It would take ",num2str(j)," years to fill a site"])
filled_sites = find(site_cap_scc_high < 0);
names_filled_sites= site_names(find(site_cap_scc_high < 0));
disp(table(filled_sites,names_filled_sites));

