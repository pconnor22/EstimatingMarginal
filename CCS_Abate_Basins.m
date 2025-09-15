%% This script creates a table to detail the basins at which CO2 is being stored

basins_num = unique(site_basin);

amnt_abate = zeros(length(basins_num),6);


sortedDPT_PV3 = table2array(macTable_scc_low(:,1));
sortedDPT45_PV3 = table2array(macTable_45Q_PV3(:,2));
sortedEmiss_PV3 = table2array(macTable_45Q_PV3(:,3));
sorted_basins_PV3 = table2array(macTable_45Q_PV3(:,7));
sortedDPT45_12_PV3 = table2array(macTable_45Q_12(:,2));



out = cell(1,numel(sortedDPT_PV3));
out1 = cell(1,numel(sortedDPT45_PV3));
out2 = cell(1,numel(sortedDPT45_12_PV3));


for i = 1:length(sortedDPT_PV3)

    if sortedDPT_PV3(i) == inf || sortedDPT_PV3(i) >= scc_low
    else
    amnt_abate(basins_num(:) == sorted_basins_PV3(i),1) = amnt_abate(basins_num(:) == sorted_basins_PV3(i),1) + sortedEmiss_PV3(i);
    out{i} = repmat(sortedDPT_PV3(i),1,sortedEmiss_PV3(i));
    end

    if sortedDPT45_PV3(i) == inf || sortedDPT45_PV3(i) > 0
    else
    amnt_abate(basins_num(:) == sorted_basins_PV3(i),2) = amnt_abate(basins_num(:) == sorted_basins_PV3(i),2) + sortedEmiss_PV3(i);
    out1{i} = repmat(sortedDPT45_PV3(i),1,sortedEmiss_PV3(i));
    end

    if sortedDPT45_12_PV3(i) == inf || sortedDPT45_12_PV3(i) > 0
    else
    amnt_abate(basins_num(:) == sorted_basins_PV3(i),3) = amnt_abate(basins_num(:) == sorted_basins_PV3(i),3) + sortedEmiss_PV3(i);
    out2{i} = repmat(sortedDPT45_12_PV3(i),1,sortedEmiss_PV3(i));
    end

end


sortedDPT_PV2 = table2array(macTable_scc_high(:,1));
% sortedDPT45_PV2 = table2array(macTable_scc_high_45Q(:,2));
sortedEmiss_PV2 = table2array(macTable_scc_high(:,2));
sorted_basins_PV2 = table2array(macTable_scc_high(:,6));

out3 = cell(1,numel(sortedDPT_PV2));
% out4 = cell(1,numel(sortedDPT45_PV2));
% out5 = cell(1,numel(sortedDPT45_PV2));


for i = 1:length(sortedDPT_PV2)

    if sortedDPT_PV2(i) == inf || sortedDPT_PV2(i) >= scc_high
    else
    amnt_abate(basins_num(:) == sorted_basins_PV2(i),4) = amnt_abate(basins_num(:) == sorted_basins_PV2(i),4) + sortedEmiss_PV2(i);
    out3{i} = repmat(sortedDPT_PV2(i),1,sortedEmiss_PV2(i));
    end

    % if sortedDPT45_PV2(i) == inf || sortedDPT45_PV2(i) > 0
    % else
    % amnt_abate(basins_num(:) == sorted_basins_PV2(i),5) = amnt_abate(basins_num(:) == sorted_basins_PV2(i),5) + sortedEmiss_PV2(i);
    % out4{i} = repmat(sortedDPT45_PV2(i),1,sortedEmiss_PV2(i));
    % end
    % 
    % if sortedDPT45_PV2(i) == inf || sortedDPT45_PV2(i) >= scc_high
    % else
    % amnt_abate(basins_num(:) == sorted_basins_PV2(i),6) = amnt_abate(basins_num(:) == sorted_basins_PV2(i),6) + sortedEmiss_PV2(i);
    % out5{i} = repmat(sortedDPT45_PV2(i),1,sortedEmiss_PV2(i));
    % end

end

macCurve = horzcat(out{:});
totalAbate = length(macCurve);
amnt_abate(:,1) = amnt_abate(:,1);

macCurve = horzcat(out1{:});
totalAbate = length(macCurve);
amnt_abate(:,2) = amnt_abate(:,2);

macCurve = horzcat(out2{:});
totalAbate = length(macCurve);
amnt_abate(:,3) = amnt_abate(:,3);

% macCurve = horzcat(out3{:});
% totalAbate = length(macCurve);
% amnt_abate(:,4) = amnt_abate(:,4);

% macCurve = horzcat(out4{:});
% totalAbate = length(macCurve);
% amnt_abate(:,5) = amnt_abate(:,5);
% 
% macCurve = horzcat(out5{:});
% totalAbate = length(macCurve);
% amnt_abate(:,6) = amnt_abate(:,6);


%geo_t=sortrows(table(basins_num,amnt_abate(:,1),amnt_abate(:,2),amnt_abate(:,3),amnt_abate(:,4),amnt_abate(:,5),amnt_abate(:,6)),2,'descend');
geo_t=sortrows(table(basins_num,amnt_abate(:,1),amnt_abate(:,2),amnt_abate(:,3),amnt_abate(:,4)),2,'descend');
geo_t(sum(table2array(geo_t(:,2:end)),2) == 0,:) = [];
% geo_t.Properties.VariableNames = ["Baisins","scc_low","45Q_PV3","45Q_12","scc_high","45Q_PV2","scc_high+45Q"];
geo_t.Properties.VariableNames = ["Baisins",num2str(scc_low),"45Q_PV3","45Q_12",num2str(scc_high)];
disp(geo_t)
writetable(geo_t,'CCS_Outputs/CCS_Tables/Geo_amnt_abate.csv')
