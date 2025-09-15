%% This script takes the outputted cost tables and assembles them into Marginal Abatement Cost Curves
% %Curves include SCC values & Abatement Levels
% close all

% dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022/total_dpt_PV3.csv');
dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/2022_SCR/total_dpt_PV3_lr.csv');
dpt_matrix = dpt_matrix(R_old,:);

%Minimizing
[dpt_array,loc] = min(dpt_matrix,[],2);

%Converting Emission data into a row vector
plant_CO2_post_capture_array=ceil(plant_CO2_post_capture/MMmt);

sites_array = site_names(loc);
basins = site_basin(loc);
st = site_state_loc(loc);


%turning the 4 row vectors into a table and then sorting
macTable =sortrows(table(dpt_array,plant_CO2_post_capture_array,plant_names,plant_fuel,sites_array,basins,st,plant_ORISPL,(1:length(dpt_array))'));

sortedDPT = table2array(macTable(:,1));
sortedEmiss = table2array(macTable(:,2));

%creating a cell array to store the $/tonne vectors for the total
%tonnes for each source
out = cell(1,numel(sortedDPT));

%looping through each source sink pair
for i = 1:length(sortedDPT)
    %excluding $/tonne costs that are infinite or egregious (> 1e8) due to
    %their distortionary effect on the curve
    if sortedDPT(i) == inf || sortedDPT(i) > 5000
    else
    %filling the specific cell of the cell array with a vector that is the
    %length of total tonnes from a source and filled  with marginal
    %cost($/tonne) of each tonne stored
    out{i} = repmat(sortedDPT(i),1,sortedEmiss(i));
    end
end

%convertingg the cell arrray into one long double row vector
macCurve = horzcat(out{:});

%plotting the curve
plot(macCurve,'LineWidth',5)
hold on
xlabel('Abatement (MMmt CO_2)','Interpreter','tex')
ylabel('Marginal Abatement Cost ($/mtCO_2)','Interpreter','tex')
ax = gca;
ax.XAxis.TickLabelFormat = '%,g';
ax.YAxis.TickLabelFormat = '%,g';
%legend('MAC','Location','northwest')
%title('Marginal Abatement Curve')
hold off
saveas(gcf,"CCS_Outputs/CCS_MAC/General_MAC.png")
writetable(macTable,"CCS_Outputs/CCS_MAC/General_MAC.xlsx")





























% %% This script takes the outputted cost tables and assembles them into Marginal Abatement Cost Curves
% %Curves include SCC values & Abatement Levels
% close all
% 
% dpt_matrix = readmatrix('CCS_Outputs/CCS_Tables/total_dpt_nom.csv');
% % dpt_matrix = readmatrix('CCS_Outputs/total_dpt_PV3.csv');
% % dpt_matrix = readmatrix('CCS_Outputs/total_dpt_PV2.csv');
% dpt_matrix = dpt_matrix(R_old,:);
% 
% %Converting Matrix into a row vector
% dpt_array = (dpt_matrix')';
% dpt_array = dpt_array(:)';
% 
% %Converting Emission data into a row vector
% plant_CO2_post_capture_matrix=(ceil(plant_CO2_post_capture/MMmt).*ones(size(dpt_matrix,1),size(dpt_matrix,2)))';
% plant_CO2_post_capture_array = plant_CO2_post_capture_matrix(:)';
% 
% %Converting Name matricies into a row vector
% plantName_matrix = strings(size(dpt_matrix,1),size(dpt_matrix,2));
% plantName_matrix = plantName_matrix + plant_names;
% plantName_array = plantName_matrix';
% plantName_array = plantName_array(:)';
% 
% site_matrix = strings(size(dpt_matrix,1),size(dpt_matrix,2));
% site_matrix = site_matrix + site_names';
% site_array = site_matrix';
% site_array = site_array(:)';
% 
% %turning the 4 row vectors into a table and then sorting
% macTable = sortrows(table(dpt_array',plant_CO2_post_capture_array',plantName_array',site_array'));
% sortedDPT = table2array(macTable(:,1));
% sortedEmiss = table2array(macTable(:,2));
% 
% %creating a cell array to store the $/tonne vectors for the total
% %tonnes for each source
% out = cell(1,numel(sortedDPT));
% 
% %looping through each source sink pair
% for i = 1:length(sortedDPT)
%     %excluding $/tonne costs that are infinite or egregious (> 1e8) due to
%     %their distortionary effect on the curve
%     if sortedDPT(i) == inf || sortedDPT(i) > 1e4
%     else
%     %filling the specific cell of the cell array with a vector that is the
%     %length of total tonnes from a source and filled  with marginal
%     %cost($/tonne) of each tonne stored
%     out{i} = repmat(sortedDPT(i),1,sortedEmiss(i));
%     end
% end
% 
% %convertingg the cell arrray into one long double row vector
% macCurve = horzcat(out{:});
% 
% %plotting the curve
% plot(macCurve,'LineWidth',5)
% xlabel('Abatement (Megatonnes CO2)')
% ylabel('Marginal Abatement Cost ($/Tonne)')
% %title('Marginal Abatement Curve')
% 
% saveas(gcf,"CCS_Outputs/CCS_MAC/General_MAC_w_Every_Pair.png")
% 
% 
% 
% 
% 
