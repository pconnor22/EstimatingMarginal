function  [storageCosts_per_cost,storageCosts,lifeTimeTonnes,dpt]=storageModule(PlumeUnArmult,nomAnnualTonne,siteNum,escalation_rate,discount_rates)
format("bank")


load("Storage_Data")
operations(1) = 12;
operations(3) = 18;
pisc_closure(2) = 19;
pisc_closure(3) = 68;

%% Site Specific Information
site = GeolDB(siteNum,:); %Desired site information
stateAbr = string(table2cell(site(1,4)));  %Storage site state abrevation
AForm = table2array(site(1,14)); %Total Surface area of injection formation
htot_ft = table2array(site(1,18)); %Thickness of injection formation
npor = table2array(site(1,21))/100; %Porosity
tmp_top = str2double(table2array(site(1,20))); %Temperature of injection formation at top of formation
Ltop = table2array(site(1,17)); %Depth to top of injection formation
kh = table2array(site(1,24)); %Horizontal permeability
litDepEnv = lower(string(table2cell(site(1,11))))+"-"+lower(string(table2cell(site(1,12))))+"-"+"reg_dip"; %Lithology-Depositional Environment


%% Calculating the CO2 Plume Area
%ACO2plz = CO2 Plume Area if plume is determined by nominal CO2 injection rate (InjProjCon = 0)
%ACO2plUn = CO2 Plume Uncertainty area
%lifeTimeTonnes = Total mass of CO2 storied in actual CO2 plume
%mCO2yr = actual mass rate of CO2 injection tonne/yr
[ACO2plz,ACO2plUn,lifeTimeTonnes,mCO2yr]=CO2PlumeArea(Aprojmaxnom,AForm, AmaxnomCon,PlumeUnArmult,htot_ft,npor,Ltop,tmp_top,...
    litDepEnv,Ang3DSeis,PlumePrFAORmult,nomAnnualTonne,operations(1),tableCoeff);

%% Calculating the Number of Vertical Wells - Active Injection Wells & Total Wells - using Law and Bachu Method
Nwell_actv = activeInjWell(lifeTimeTonnes,operations(1),maxCO2mult,htot_ft,Ltop, kh, tmp_top,qwell_mech_day,grad_Plith,fact_Pfrac);
Nwell_fin = vertInjWell(Nwell_actv);

%% Calculating Site Specific Costs
cost_values = spc(cost_values,stateAbr,site_selection_characterization,operations,permitting_construction,pisc_closure,drillCosts,project_contingency,process_contingency,rathole,addDepthStrat,Ltop,htot_ft,...
    groundwater_depth,Ppumpin,Ppumpout,maxCO2mult,ACO2plz,lifeTimeTonnes,Nwell_fin,Nwell_actv,ACO2plUn,dens_wells_CA,PlumePrFAORmult,apipewell,mCO2yr,PlumeUnArmult);

%% Plume and Well Schedule
[plume_area_other,monitoring_wells_reservoir,monitoring_wells_dual,monitoring_wells_above,monitoring_ground_water,monitoring_vadose] = plume_well_schedule(site_selection_characterization,permitting_construction,operations, pisc_closure,...
    PlumeUnArmult,ACO2plz,PlumePrFAORmult,Ang3DSeis,Ang2DSeis,frac_3D,frac_2D,ACO2plUn,Ltop,htot_ft,num_strat_test_wells,...
    strat_2_inj,Nwell_fin,well_spacing);

%% Cashflow
cf=cashflows(pisc_closure,cost_values,cost_categorized,wc_depend,Nwell_fin,plume_area_other,monitoring_wells_reservoir,monitoring_wells_dual,monitoring_wells_above,monitoring_ground_water,monitoring_vadose...
    ,pl_depend,t_depend);

diff = start_year - base_year;
esc_sched = zeros(1,pisc_closure(3));

i=1:pisc_closure(3);
esc_sched(i) = (1+escalation_rate).^(diff+i-1);
disc_sched_3(i) = 1./(1+discount_rates(1,1)).^(i-1);
disc_sched_2(i) = 1./(1+discount_rates(1,2)).^(i-1);

final_cf = esc_sched.*sum(cf);
dpt_nom=sum(final_cf)/lifeTimeTonnes;
dpt_PV_3 = sum(disc_sched_3.*final_cf)/lifeTimeTonnes;
dpt_PV_2 =  sum(disc_sched_2.*final_cf)/lifeTimeTonnes;

storageCosts_per_cost = sum(cf,2);
storageCosts = final_cf;
dpt = [dpt_nom;dpt_PV_3;dpt_PV_2];

end

%% Cashflows
function cf=cashflows(pisc_closure,cost_values,cost_categorized,wc_depend,Nwell_fin,plume_area_other,monitoring_wells_reservoir,monitoring_wells_dual,monitoring_wells_above,monitoring_ground_water,monitoring_vadose,pl_depend,t_depend)

i = 1:pisc_closure(3);
cf = zeros(length(cost_values),length(pisc_closure(3)));

for j = 1:length(cost_values)

    if wc_depend(j) == "yes - total inj wells (reg & convrt)" || wc_depend(j) ==  "yes - total inject wells (reg & convert)"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1) * Nwell_fin;

    elseif wc_depend(j) == "yes - total deep mon wells (reg & convrt; in res, abv seal & dual)"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1) * (max(monitoring_wells_reservoir(13,:))+ max(monitoring_wells_dual(13,:)) + max(monitoring_wells_above(13,:)));

    elseif wc_depend(j) == "yes - total wat disp wells"
        % This needs to be modified if you ever turn water on, needs mult
        % by number of water wells
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1) .* plume_area_other(24,i);

    elseif wc_depend(j) == "yes - total wat prod wells"
        % This needs to be modified if you ever turn water on, needs mult
        % by number of water wells
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1) .* plume_area_other(22,i);

    elseif wc_depend(j) == "yes - total in res mon wells (reg & convrt)" || wc_depend(j) == "yes - total in res mon wells (reg & convert)"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_reservoir(13,i);

    elseif wc_depend(j) == "yes - total dual compl mon wells (reg & convrt)" || wc_depend(j) == "yes - total dual compl mon wells (reg & convert)"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_dual(13,i);

    elseif wc_depend(j) == "yes - total ab seal mon wells (reg & convrt)" || wc_depend(j) == "yes - total abv seal mon wells (reg & convert)" || wc_depend(j) == "yes - total ab seal mon wells (reg & convert)" || wc_depend(j) == "yes - total abv seal mon wells (reg & convrt)"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_above(13,i);

    elseif wc_depend(j) == "yes - total GW mon wells (reg & convrt)" || wc_depend(j) == "yes - total GW mon wells (reg & convert)" || wc_depend(j) == "yes - total GW mon wells"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_ground_water(2,i);

    elseif wc_depend(j) == "yes - total vad zone wells" || wc_depend(j) =="yes - total vad zone mon wells"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_vadose(2,i);

    elseif wc_depend(j) == "yes - total reg strat wells"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*plume_area_other(11,i);


    elseif wc_depend(j) == "yes - new all strat wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*plume_area_other(8,i);

    elseif wc_depend(j) == "yes - new reg inj wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*plume_area_other(17,i);

    elseif wc_depend(j) == "yes - new reg in res mon wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_reservoir(12,i);

    elseif wc_depend(j) == "yes - new reg above seal mon wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_above(12,i);

    elseif wc_depend(j) == "yes - new reg dual compl mon wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_dual(12,i);

    elseif wc_depend(j) == "yes - new GW mon wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_ground_water(1,i);

    elseif wc_depend(j) == "yes - new vad zone mon wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_vadose(1,i);

    elseif wc_depend(j) == "yes - new wat prod wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*plume_area_other(21,i);

    elseif wc_depend(j) == "yes - new wat disp wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*plume_area_other(23,i);

    elseif wc_depend(j) == "yes - total wat prod & wat disp wells"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*(plume_area_other(22,i)+plume_area_other(24,i));


    elseif wc_depend(j) == "yes - convert strat to inj wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*plume_area_other(15,i);

    elseif wc_depend(j) == "yes - convert strat to in res mon wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_reservoir(11,i);

    elseif wc_depend(j) == "yes - convert strat to above seal mon wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_above(11,i);

    elseif wc_depend(j) == "yes - convert strat to dual compl mon wells, additional"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1).*monitoring_wells_dual(11,i);


    elseif pl_depend(j) == "yes -3D seismic"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1) .* plume_area_other(5,i);

    elseif pl_depend(j) == "yes -3D AoR"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1) .* plume_area_other(2,i);

    elseif pl_depend(j) == "yes -2D seismic"
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1) .* plume_area_other(7,i);


    else
        cf(j,i) = ((cost_values(j,2) <= i & i <= cost_values(j,3)) & (strcmp(cost_categorized(j), 'Monitoring')...
            & (i == cost_values(j,3) | mod((i - (cost_values(j,2) + cost_values(j,4) - 1)), cost_values(j,4)) == 0)...
            | mod((i - cost_values(j,2)), cost_values(j,4)) == 0)) * cost_values(j,1);

    end


end



end


%% Number of Active Injection wells & Helper Functions
function Nwell_actv = activeInjWell(lifeTimeTonnes,yrsOp,maxCO2mult,htot_ft,Ltop, kh, tmp_top,qwell_mech_day,grad_Plith,fact_Pfrac)


mCO2maxdy = (lifeTimeTonnes/yrsOp)*maxCO2mult/365;
grad_Phyd = 0.464; %typical value for brine, ambient hydrostatic pressure gradient
Pamb_top = grad_Phyd * Ltop;
Pamb_mid = Pamb_top + 0.5 * htot_ft * grad_Phyd;%Psia
Pamb_mid = 0.0068948 * Pamb_mid;%converting to MPaPamb_mid = Pamb_top + 0.5 * htot_ft * grad_Phyd
coef_LB =  0.0208;
fact_kv_h	= 0.3;
fact_Pinj	= 0.90;


Pfrac_top = fact_Pfrac * grad_Plith * Ltop;


Phydmax_top = fact_Pinj * Pfrac_top;
Phydmax_mid = Phydmax_top + 0.5*htot_ft*grad_Phyd;
Phydmax_mid = 0.0068948  * Phydmax_mid;


kv = kh * fact_kv_h;
kabs = sqrt(kh * kv);

grad_tmp	= 1.37;
tmp_mid = tmp_top + 0.5*htot_ft*(0.01)*grad_tmp;
tmp_mid = 273.15 + (5/9.0) * (tmp_mid-32);


visCO2amb_mid_uPa =visPRCO2(denPRCO2(0,Pamb_mid,tmp_mid),tmp_mid);
visCO2amb_mid_cp = visCO2amb_mid_uPa * 0.001;

visCO2inj_mid_uPa =visPRCO2(denPRCO2(0,Phydmax_mid,tmp_mid),tmp_mid);
visCO2inj_mid_cp = visCO2inj_mid_uPa * 0.001;
visCO2av_mid = 0.5 * (visCO2amb_mid_cp + visCO2inj_mid_cp);

fact_hn_tot = 1.00;
hnet = fact_hn_tot * htot_ft;
hnet = hnet / (3.281);

qinj_LB = coef_LB * kabs / visCO2av_mid;
qwell_form_LB = qinj_LB * hnet * (Phydmax_mid - Pamb_mid);
qwell_LB = min(qwell_form_LB, qwell_mech_day);



Nwell_actv = ceil(mCO2maxdy / qwell_LB);

end

%Calculates Viscosity of CO2
function vis = visPRCO2(den,tmp)

tmpcs = 251.196;
az(1) = 0.235156;
az(2) = -0.491266;
az(3) = 0.05211155;
az(4) = 0.05347906;
az(5) = -0.01537102;

CSR = 0;
for j = 1:5
    CSR = CSR + az(j) * (log(tmp / tmpcs)) ^ (j - 1);
end

visz = 1.00697 * sqrt(tmp) / exp(CSR);


%Calculation of excess viscosity term that corrects
% for non-zero density
d11 = 0.004071119;
d21 = 7.198037E-05;
d64 = 2.411697E-17;
d81 = 2.971072E-23;
d82 = -1.627888E-23;

visexc = d11 * den + d21 * den ^ 2 + d64 * den ^ 6 * (tmp / tmpcs) ^ -3 + d81 * den ^ 8 + d82 * den ^ 8 * (tmp / tmpcs) ^ -1;

%Calculate viscosity
vis = visz + visexc;
end

function Nwell_fin_LB = vertInjWell(Nwell_actv_LB)

fact_well_add = .10;

Nwell_fin_a = ceil(Nwell_actv_LB * (1+fact_well_add));
Nwell_min = 4;

if Nwell_fin_a <=  Nwell_min
    Nwell_fin_LB = Nwell_min;
else
    Nwell_fin_LB = Nwell_fin_a;
end

end


%% Actual Annual Rate of CO2 Injection
function [mCO2yr,mCO2tot] = annualMassCO2(nomAnnualTonne,Injdur,ACO2pl,npor,Ltop,tmp_top,litDepEnv,htot_ft,tableCoeff)

volinit = 0;
grad_Phyd = 0.464; %typical value for brine, ambient hydrostatic pressure gradient
Pamb_top = grad_Phyd * Ltop;
Pamb_mid = Pamb_top + 0.5 * htot_ft * grad_Phyd;%Psia
Pamb = 0.0068948 * Pamb_mid;%converting to MPa
pres = Pamb;

grad_tmp = 1.37;%typical value for temperature gradient, degF/ft
tmp_surf = 59; %typical value for surface tmperature, degF
if isnan(tmp_top)
    tmp_top = tmp_surf + grad_tmp * Ltop * 0.01;
end
tmp_mid = tmp_top + 0.5*htot_ft*0.01*grad_tmp;%temperataure at midpoint injection formation defF
tmp= 273.15 + (5.0/9) * (tmp_mid-32);%converting to Kelvin


lith =table2array(tableCoeff(:,1));
coeff =table2array(tableCoeff(:,2));
E = coeff(litDepEnv == lith);

denCO2=denPRCO2(volinit, pres, tmp);

mCO2nomtot = nomAnnualTonne * Injdur;
ACO2plnom = (mCO2nomtot * 1000)/ (htot_ft*0.3048 * npor * denCO2 * E);
ACO2plnom = ACO2plnom / 2589988.110336;

mCO2tot = mCO2nomtot * ACO2pl / ACO2plnom;
mCO2yr = mCO2tot / Injdur;


end


%% CO2 Plume Area Function & Helper Functions
%AForm varies for formation
%Aprojmaxnom is always defaulted to 100 at this point key inputs
%AmaxnomCon should be 0, default in model sets CO2 size
%PlumeUnArmult, plume uncertainty area, defaulted to 1.75
%Htot height of injection, specific to formation
%npor best estimate porosity, specific to formation
%Ltop Depth to top of injection formation, specific to formation mi2
%tmp_top variable may be availble specific to formation, otherwise
%calculated, selected temperaature at top of injection formation degF\
%LithDepEnv tells us the ltihology depositional environment to find the
%storage coefficent, it is specific to the formation, all have reg_dip to
%the assumption
%ang3Dseis default parameter set at 45 degrees
%plUnMult plume uncertaintity multiplier set at a default value of 1.75
%OUTPUTS Max CO2 Plume AREA mi2
function [ACO2plz,ACO2plUn,lifeTimeTonnes,mCO2yr]=CO2PlumeArea(Aprojmaxnom,AForm, AmaxnomCon,PlumeUnArmult,htot_ft,npor,Ltop,tmp_top,litDepEnv,Ang3DSeis,PlumePrFAORmult,nomAnnualTonne,yrOp,tableCoeff)

volinit = 0;
grad_Phyd = 0.464; %typical value for brine, ambient hydrostatic pressure gradient
Pamb_top = grad_Phyd * Ltop;
Pamb_mid = Pamb_top + 0.5 * htot_ft * grad_Phyd;%Psia
Pamb = 0.0068948 * Pamb_mid;%converting to MPa
pres = Pamb;

grad_tmp = 1.37;%typical value for temperature gradient, degF/ft
tmp_surf = 59; %typical value for surface tmperature, degF
if isnan(tmp_top)
    tmp_top = tmp_surf + grad_tmp * Ltop * 0.01;
end
tmp_mid = tmp_top + 0.5*htot_ft*0.01*grad_tmp;%temperataure at midpoint injection formation defF
tmp= 273.15 + (5.0/9) * (tmp_mid-32);%converting to Kelvin

lith =table2array(tableCoeff(:,1));
coeff =table2array(tableCoeff(:,2));
E = coeff(litDepEnv == lith);

denCO2=denPRCO2(volinit, pres, tmp);
mCO2nomtot = nomAnnualTonne * yrOp;
ACO2plnom = mCO2nomtot * 1000 / (htot_ft*0.3048 * npor * denCO2 * E);
ACO2plnom = ACO2plnom / 2589988.110336;


PForm_struc = 0.9750;%set to values due to assuming reg dip structure
PAvail_struc = 0.400;

AFStructmax = AForm * PForm_struc * PAvail_struc;
Aprojmax = min([Aprojmaxnom, AFStructmax]);

Lbot = (Ltop + htot_ft)/5280; %feet converted to miles,depth to bottom of formation
ACO2plmax1 = Aprojmax / PlumeUnArmult;
ACO2plmax2 = Aprojmax / (PlumeUnArmult * PlumePrFAORmult);
ACO2plmax3 = pi * (sqrt(Aprojmax / pi) - Lbot - tand(Ang3DSeis))^2 / PlumeUnArmult;

if( AmaxnomCon == 1)
    if(ACO2plnom > ACO2plmax3)
        ACO2plz = ACO2plmax3;
    else
        ACO2plz = ACO2plnom;
    end
elseif AmaxnomCon == 2
    if ACO2plnom >min([ACO2plmax2, ACO2plmax3])
        ACO2plz = min([ACO2plmax2, ACO2plmax3]);
    else
        ACO2plz = ACO2plnom;
    end
else
    if(ACO2plnom > ACO2plmax1)
        ACO2plz = ACO2plmax1;
    else
        ACO2plz = ACO2plnom;
    end
end

[mCO2yr,lifeTimeTonnes] = annualMassCO2(nomAnnualTonne,yrOp,ACO2plz,npor,Ltop,tmp_top,litDepEnv,htot_ft,tableCoeff);

ACO2plUn = ACO2plz * PlumeUnArmult;

end

%Calculates the density of CO2
function den=denPRCO2(volinit, pres, tmp)

MW = 44.0095; %molecular weight of CO2 (g/mol)

%Calculate molar volume in m3/mol, then density

den = MW * 0.001 / volPRCO2(volinit, pres, tmp);

end

function volly=volPRCO2(volinit, pres, tmp)


% This function calculates the molar volume
% resulting from the Peng-Robinson equation of state for CO2.

%volinit is an inital guess of the molar volume
% pres is the absolute pressure (MPa)
% tmp is the temperature (deg K)

% The function returns the molar volume in m3/mol.


Rgas = 8.314E-06;    %universal gas constant (m3-MPa/degK-mol)
tol = 1E-05;   %desired accuracy of result

%Properties of CO2
tmpcrit = 304.1282;   % critical temperature (deg K)
prescrit = 7.3773; % critical pressure (MPa)

b1 = 0.0778 * Rgas * tmpcrit / prescrit;


%Determine inital guess
if (volinit <= 0)
    vol1a = Rgas * tmp / pres;
    if (vol1a < b1)
        vol1a = b1 * 1.1;
    end
else
    vol1a = volinit;
end


%Determine molar volume using Newton-Raphson method
%The inner loop performs Newton Raphson method
%Limit inner loop to 1000 iterations
%The outer loop checks for volume being positive.
%If volume is negative, intial volume was too low.
%Reset intial volume to increase value (1.05 times previous
%initial volume)
%Rerun and check until it converges
%Limit outer loop to 100 iterations

for i = 1:100
    if(i == 1)

    else
        if (vol2 > 0)
            break;
        else
            vol1a = vol1a * 1.05;
        end
    end
    vol1 = vol1a;
    for j = 1:1000
        vol2 = vol1 - fvolPRCO2(vol1, pres, tmp) / d1fvolPRCO2(vol1, pres, tmp);
        reldif = abs(vol1 - vol2) / abs(vol2);
        if (reldif <= tol)
            break;
        else
            vol1 = vol2;
        end
    end
end

volly = vol2;

end


function fvolly=fvolPRCO2(vol, pres, tmp)

%This function calculates the polynomial f(vol,pres,tmp)
% resulting from the Peng-Robinson equation of state for CO2.
% When f(vol,pres,tmp)=0, then vol, pres and tmp are a solution
% of the the Peng-Robinson equation of state for CO2.

% vol is the molar volume (m3/mol)
% pres is the absolute pressure (MPa)
% tmp is the temperature (deg K)

% The function returns f(vol, pres,tmp).



Rgas = 8.314E-06;    % universal gas constant (m3-MPa/degK-mol)

% Properties of CO2
tmpcrit = 304.1282;   % critical temperature (deg K)
prescrit = 7.3773;  %critical pressure (MPa)
waf = 0.22394;  %accentric factor

m1 = 0.37464 + 1.54226 * waf - 0.26992 * waf * waf;
a1 = 0.45724 * (Rgas * tmpcrit) ^ 2 * (1 + m1 * (1 - sqrt(tmp / tmpcrit))) ^ 2 / prescrit;
b1 = 0.0778 * Rgas * tmpcrit / prescrit;

fvolly = pres * vol ^ 3 + (pres * b1 - Rgas * tmp) * vol ^ 2 + (a1 - 3 * pres * b1 ^ 2 - 2 * Rgas * tmp * b1) * vol + pres * b1 ^ 3 + Rgas * tmp * b1 ^ 2 - a1 * b1;

end


function d1volly=d1fvolPRCO2(vol, pres, tmp)


%This function calculates the first derivative with respect to vol
% of the polynomial f(vol,pres,tmp)
% resulting from the Peng-Robinson equation of state for CO2.
% This function is used in the Newton-Raphson method
% to find a root of the equation f(vol,pres,tmp)=0,
% given values for pres and tmp.

% vol is the molar volume (m3/mol)
% pres is the absolute pressure (MPa)
% tmp is the temperature (deg K)
%
%The function returns d f(vol, pres,tmp) / d vol.



Rgas = 8.314E-06; % universal gas constant (m3-MPa/degK-mol)

% Properties of CO2
tmpcrit = 304.1282; % critical temperature (deg K)
prescrit = 7.3773; % critical pressure (MPa)
waf = 0.22394;  %accentric factor

m1 = 0.37464 + 1.54226 * waf - 0.26992 * waf * waf;
a1 = 0.45724 * (Rgas * tmpcrit) ^ 2 * (1 + m1 * (1 - sqrt(tmp / tmpcrit))) ^ 2 / prescrit;
b1 = 0.0778 * Rgas * tmpcrit / prescrit;

d1volly = 3 * pres * vol ^ 2 + 2 * (pres * b1 - Rgas * tmp) * vol + (a1 - 3 * pres * b1 ^ 2 - 2 * Rgas * tmp * b1);

end


function cost_values = spc(cost_values,stateAbr,site_selection_characterization,operations,permitting_construction,pisc_closure,drillCosts,project_contingency,process_contingency,rathole,addDepthStrat,Ltop,htot_ft,...
    groundwater_depth,Ppumpin,Ppumpout,maxCO2mult,ACO2plz,lifeTimeTonnes,Ninjwell,Ninjwellactiv,ACO2plUn,dens_wells_CA,PlumePrFAORmult,apipewell,mCO2yr,PlumeUnArmult)


%% Well-Drilling Specific Costs

%Well timelines
strat=[site_selection_characterization(2),site_selection_characterization(2),1];
inj = [operations(2),operations(3),1];
monitoring = [operations(2),operations(3),1];


states = string(table2cell(drillCosts(:,1)));

%4.2 Drilling Costs
cost_values(306,2) = strat(1);
cost_values(306,3) = strat(2);
cost_values(306,4) = 1;

cost = [table2array(drillCosts(states==stateAbr,2)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(306,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(307,2) = permitting_construction(3);
cost_values(307,3) = permitting_construction(3);
cost_values(307,4) = 1;

cost = [table2array(drillCosts(states==stateAbr,3)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(307,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(308,2) = permitting_construction(3);
cost_values(308,3) = permitting_construction(3);
cost_values(308,4) = 1;

if table2array(drillCosts(states==stateAbr,3)) - table2array(drillCosts(states==stateAbr,2)) > 0
    cost = [table2array(drillCosts(states==stateAbr,3)) - table2array(drillCosts(states==stateAbr,2)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
    cost_values(308,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
else
    cost = [0,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
    cost_values(308,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
end

cost_values(309,2) = monitoring(1);
cost_values(309,3) = monitoring(2);
cost_values(309,4) = 1;

cost = [table2array(drillCosts(states==stateAbr,4)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(309,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


cost_values(310,2) = monitoring(1);
cost_values(310,3) = monitoring(2);
cost_values(310,4) = 1;

if table2array(drillCosts(states==stateAbr,4)) - table2array(drillCosts(states==stateAbr,2)) > 0
    cost = [table2array(drillCosts(states==stateAbr,4)) - table2array(drillCosts(states==stateAbr,2)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
    cost_values(310,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
else
    cost = [0,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
    cost_values(310,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
end

cost_values(311,2) = monitoring(1);
cost_values(311,3) = monitoring(2);
cost_values(311,4) = 1;

cost = [table2array(drillCosts(states==stateAbr,5)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(311,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


cost_values(312,2) = monitoring(1);
cost_values(312,3) = monitoring(2);
cost_values(312,4) = 1;

if table2array(drillCosts(states==stateAbr,5)) - table2array(drillCosts(states==stateAbr,2)) > 0
    cost = [table2array(drillCosts(states==stateAbr,5)) - table2array(drillCosts(states==stateAbr,2)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
    cost_values(312,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
else
    cost = [0,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
    cost_values(312,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
end

cost_values(313,2) = monitoring(1);
cost_values(313,3) = monitoring(2);
cost_values(313,4) = 1;


cost = [table2array(drillCosts(states==stateAbr,4)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(313,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(314,2) = monitoring(1);
cost_values(314,3) = monitoring(2);
cost_values(314,4) = 1;

if table2array(drillCosts(states==stateAbr,4)) - table2array(drillCosts(states==stateAbr,2)) > 0
    cost = [table2array(drillCosts(states==stateAbr,4)) - table2array(drillCosts(states==stateAbr,2)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
    cost_values(310,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
else
    cost = [0,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
    cost_values(314,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
end

cost_values(315,2) = monitoring(1);
cost_values(315,3) = monitoring(2);
cost_values(315,4) = 1;

cost = [table2array(drillCosts(states==stateAbr,6)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(315,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(316,2) = monitoring(1);
cost_values(316,3) = monitoring(2);
cost_values(316,4) = 1;

cost = [table2array(drillCosts(states==stateAbr,7)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(316,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(317,2) = monitoring(1);
cost_values(317,3) = monitoring(2);
cost_values(317,4) = 1;

cost = [table2array(drillCosts(states==stateAbr,8)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(317,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(318,2) = monitoring(1);
cost_values(318,3) = monitoring(2);
cost_values(318,4) = 1;

cost = [table2array(drillCosts(states==stateAbr,9)),1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(318,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%319 - 331 =0

%4.3 Wireline (Geophysical) Logging
Lbot = (Ltop + htot_ft);

strat_depth = Lbot + rathole + addDepthStrat;
injection_depth = Ltop + htot_ft + rathole;
inReservoir_depth = Ltop + htot_ft + rathole;
aboveSeal_depth = Ltop - 200 + rathole;
dualCompletion_depth = inReservoir_depth;

waterProduction_depth = injection_depth;
waterDisposal_depth = Ltop;


cost_values(332,2) = strat(1);
cost_values(332,3) = strat(2);
cost_values(332,4) = 1;

cost = [strat_depth*2.75,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(332,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(333,2) = permitting_construction(3);
cost_values(333,3) = permitting_construction(3);
cost_values(333,4) = 1;


cost = [injection_depth*7.15+2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(333,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(334,2) = permitting_construction(3);
cost_values(334,3) = permitting_construction(3);
cost_values(334,4) = 1;

cost = [2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(334,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(335,2) = monitoring(1);
cost_values(335,3) = monitoring(2);
cost_values(335,4) = 1;

cost = [inReservoir_depth*7.15+2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(335,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(336,2) = monitoring(1);
cost_values(336,3) = monitoring(2);
cost_values(336,4) = 1;

cost = [2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(336,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


cost_values(337,2) = monitoring(1);
cost_values(337,3) = monitoring(2);
cost_values(337,4) = 1;

cost = [aboveSeal_depth*5.65+2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(337,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(338,2) = monitoring(1);
cost_values(338,3) = monitoring(2);
cost_values(338,4) = 1;

cost = [2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(338,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(339,2) = monitoring(1);
cost_values(339,3) = monitoring(2);
cost_values(339,4) = 1;

cost = [dualCompletion_depth*7.15+2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(339,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(340,2) = monitoring(1);
cost_values(340,3) = monitoring(2);
cost_values(340,4) = 1;

cost = [2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(340,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(341,2) = monitoring(1);
cost_values(341,3) = monitoring(2);
cost_values(341,4) = 1;

cost = [groundwater_depth*0.75,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(341,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(342,2) = monitoring(1);
cost_values(342,3) = monitoring(2);
cost_values(342,4) = 1;

cost = [waterProduction_depth*5.65+2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(342,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(343,2) = monitoring(1);
cost_values(343,3) = monitoring(2);
cost_values(343,4) = 1;

cost = [waterDisposal_depth*5.65+2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(343,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%4.4 Core Recovery 352
cost_values(344,2) = strat(1);
cost_values(344,3) = strat(2);
cost_values(344,4) = 1;

cost_values(345,2) = permitting_construction(3);
cost_values(345,3) = permitting_construction(3);
cost_values(345,4) = 1;

cost_values(346:352,2) = monitoring(1);
cost_values(346:352,3) = monitoring(2);
cost_values(346:352,4) = 1;

cost = [18100,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(344:348,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
cost_values(351:352,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%4.5 Fluid Recovery
cost_values(353,2) = strat(1);
cost_values(353,3) = strat(2);
cost_values(353,4) = 1;

cost_values(354,2) = permitting_construction(3);
cost_values(354,3) = permitting_construction(3);
cost_values(354,4) = 1;

cost_values(355:360,2) = monitoring(1);
cost_values(355:360,3) = monitoring(2);
cost_values(355:360,4) = 1;

cost = [2207,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(353:357,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
cost = [1207,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(358,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%359 - 360 = 0

%3.6 Well Tests
cost_values(361,2) = strat(1);
cost_values(361,3) = strat(2);
cost_values(361,4) = 1;

cost_values(362:363,2) = permitting_construction(3);
cost_values(362:363,3) = permitting_construction(3);
cost_values(362:363,4) = 1;

cost_values(364:372,2) = monitoring(1);
cost_values(364:372,3) = monitoring(2);
cost_values(364:372,4) = 1;

cost = [2070,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(361:2:365,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
cost_values(369,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [4140,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(368,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);
cost_values(362:2:364,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%4.7 Well Seismic
cost_values(373,2) = strat(1);
cost_values(373,3) = strat(2);
cost_values(373,4) = 1;

cost_values(374,2) = permitting_construction(3);
cost_values(374,3) = permitting_construction(3);
cost_values(374,4) = 1;

cost_values(375:379,2) = monitoring(1);
cost_values(375:379,3) = monitoring(2);
cost_values(375:379,4) = 1;

cost = [300000,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(373:377,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%378-379 = 0

%4.8 Analysis
%380 - 388 = 0

%4.9 Completion
diameter_casing_in = 7.50;
diameter_tubing_in = 4;

cost_values(389:390,2) = inj(1);
cost_values(389:390,3) = inj(2);
cost_values(389:390,4) = 1;

cost = [15000,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(389:390,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%391 - 399 = 0
cost_values(400:401,2) = permitting_construction(3);
cost_values(400:401,3) = permitting_construction(3);
cost_values(400:401,4) = 1;

cost = [2.7 * diameter_casing_in * injection_depth + 1.15 * diameter_tubing_in*injection_depth + 3 * injection_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(400:401,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(402:403,2) = monitoring(1);
cost_values(402:403,3) = monitoring(2);
cost_values(402:403,4) = 1;

cost = [2.7 * diameter_casing_in * inReservoir_depth + 3 * inReservoir_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(402:403,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%404 - 405 = 0

cost_values(406:407,2) = monitoring(1);
cost_values(406:407,3) = monitoring(2);
cost_values(406:407,4) = 1;

cost = [2.7 * diameter_casing_in * dualCompletion_depth + 3 * dualCompletion_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(406:407,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%408 - 409 = 0
cost_values(408:409,2) = monitoring(1);
cost_values(408:409,3) = monitoring(2);
cost_values(408:409,4) = 1;

% 4.10 Downhole Equipment for Wells
cost_values(410:411,2) = permitting_construction(3);
cost_values(410:411,3) = permitting_construction(3);
cost_values(410:411,4) = 1;

cost = [500,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(410:411,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(412:419,2) = monitoring(1);
cost_values(412:419,3) = monitoring(2);
cost_values(412:419,4) = 1;

cost = [10400,1,0,0,0,0,0,0,0,0,1,1,1+process_contingency,1+project_contingency];
cost_values(412:415,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [20800,1,0,0,0,0,0,0,0,0,1,1,1+process_contingency,1+project_contingency];
cost_values(416:417,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%4.20 Operations & Maintenance
%New Timing

% strat=[site_selection_characterization(2),site_selection_characterization(3)];
inj = [operations(2),operations(3)];
monitoring = [operations(2),operations(3)];
monitoring_pisc = [pisc_closure(2),pisc_closure(3)];

cost_values(420,2) = inj(1);
cost_values(420,3) = inj(2);
cost_values(420,4) = 1;


cost = [77500+3.1*injection_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(420,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(421,2) = monitoring(1);
cost_values(421,3) = monitoring(2);
cost_values(421,4) = 1;

cost_values(422,2) = monitoring_pisc(1);
cost_values(422,3) = monitoring_pisc(2);
cost_values(422,4) = 1;

cost = [25900+3.1*inReservoir_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(421:422,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(423,2) = monitoring(1);
cost_values(423,3) = monitoring(2);
cost_values(423,4) = 1;

cost_values(424,2) = monitoring_pisc(1);
cost_values(424,3) = monitoring_pisc(2);
cost_values(424,4) = 1;

cost = [25900+3.1*aboveSeal_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(423:424,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(425,2) = monitoring(1);
cost_values(425,3) = monitoring(2);
cost_values(425,4) = 1;

cost_values(426,2) = monitoring_pisc(1);
cost_values(426,3) = monitoring_pisc(2);
cost_values(426,4) = 1;

cost = [25900+3.1*dualCompletion_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(425:426,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(427,2) = monitoring(1);
cost_values(427,3) = monitoring(2);
cost_values(427,4) = 1;

cost_values(428,2) = monitoring_pisc(1);
cost_values(428,3) = monitoring_pisc(2);
cost_values(428,4) = 1;

cost = [2000,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(427:428,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(429,2) = monitoring(1);
cost_values(429,3) = monitoring(2);
cost_values(429,4) = 1;

cost_values(430,2) = monitoring_pisc(1);
cost_values(430,3) = monitoring_pisc(2);
cost_values(430,4) = 1;

cost = [100,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(429:430,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(431:432,2) = monitoring(1);
cost_values(431:432,3) = monitoring(2);
cost_values(431:432,4) = 1;

cost = [77500+3.1*waterProduction_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(431,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [77500+3.1*waterDisposal_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(432,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


%4.21 Mechanical Integrity Tests
cost_values(433,2) = inj(1);
cost_values(433,3) = inj(2);
cost_values(433,4) = 1;

cost = [4140+8.3*injection_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(433,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(434,2) = monitoring(1);
cost_values(434,3) = monitoring(2);
cost_values(434,4) = 5;

cost_values(435,2) = monitoring_pisc(1);
cost_values(435,3) = monitoring_pisc(2);
cost_values(435,4) = 5;


cost = [2070+4.15*inReservoir_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(434:435,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(436,2) = monitoring(1);
cost_values(436,3) = monitoring(2);
cost_values(436,4) = 5;


cost_values(437,2) = monitoring_pisc(1);
cost_values(437,3) = monitoring_pisc(2);
cost_values(437,4) = 5;

cost_values(438,2) = monitoring(1);
cost_values(438,3) = monitoring(2);
cost_values(438,4) = 5;

cost_values(439,2) = monitoring_pisc(1);
cost_values(439,3) = monitoring_pisc(2);
cost_values(439,4) = 5;

cost = [2070+4.15*dualCompletion_depth,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(438:439,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(440:441,2) = monitoring(1);
cost_values(440:441,3) = monitoring(2);
cost_values(440:441,4) = 5;

cost = [11410,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(440:441,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%4.30 Periodic Monitoring
%442-446 = 0

cost_values(447,2) = monitoring(1);
cost_values(447,3) = monitoring(2);
cost_values(447,4) = 5;

cost_values(448,2) = monitoring_pisc(1);
cost_values(448,3) = monitoring_pisc(2);
cost_values(448,4) = 5;

cost = [300000,1,0,0,0,0,0,0,0,0,1,1,1+process_contingency,1+project_contingency];
cost_values(447:448,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%4.40  - % 4.45Monitoring for Conformance
%449-452 = 0

%4.50 Plugging and Abandon

strat = [site_selection_characterization(3),site_selection_characterization(3)];
inj = [operations(3),operations(3)];
monitoring = [pisc_closure(3),pisc_closure(3)];


cost_values(453,2) = strat(1);
cost_values(453,3) = strat(2);
cost_values(453,4) = 1;

cost = [50900,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(453,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(454,2) = inj(1);
cost_values(454,3) = inj(2);
cost_values(454,4) = 1;

cost = [(115370+4.235*injection_depth),1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(454,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


cost_values(455:459,2) = monitoring(1);
cost_values(455:459,3) = monitoring(2);
cost_values(455:459,4) = 1;

cost = [91600,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(455,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [71600,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(456,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [91600,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(457,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [2000,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(458,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [700,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(459,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


cost_values(460:461,2) = inj(1);
cost_values(460:461,3) = inj(2);
cost_values(460:461,4) = 1;


cost = [112300,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(460:461,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%% Surface Equipment Costs

%Feeder Pipeline Capital
cost_values(462,2) = operations(2);
cost_values(462,3) = operations(2);
cost_values(462:476,4) = 1;


mCO2actyr = (lifeTimeTonnes/operations(1))*0.000001;

if mCO2actyr < 0.29
    Cpipfdcapvar = 360000 ;
elseif mCO2actyr < 0.676
    Cpipfdcapvar = 540000;
elseif mCO2actyr < 1.16
    Cpipfdcapvar = 720000;
elseif mCO2actyr < 5.0
    Cpipfdcapvar = 900000;
else
    Cpipfdcapvar = 1200000;
end

Lpipfd =2*sqrt(ACO2plz/pi)/2;
Cpipfdcapfix = 200000;
Cpipfdcaptot = Cpipfdcapfix + Lpipfd * Cpipfdcapvar;

cost = [Cpipfdcaptot,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(462,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Feeder Pipeline O&M
cost_values(463,2) = operations(2);
cost_values(463,3) = operations(3);

CpipfdOMvar	= 9000;
CpipfdOMtot = Lpipfd * CpipfdOMvar;

cost = [CpipfdOMtot,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(463,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Capital Costs for Pressure Boosting Pump
cost_values(464,2) = operations(2);
cost_values(464,3) = operations(2);

Cpumpcapfix	= 87000;
Cpumpcapvar	= 1400;

epump = 0.75;

tmp = 284.816666667; %kelvin
Ppumpin = 0.00689475729 * Ppumpin;
Ppumpout = 0.00689475729 * Ppumpout;
pres = 0.5 * (Ppumpin + Ppumpout);
denCO2pump=denPRCO2(0, pres, tmp);

mCO2actmaxsec = maxCO2mult * mCO2actyr * 1e9 / (365*24*3600);
Wpump = mCO2actmaxsec * (Ppumpout - Ppumpin) * 1e6 * 1e-3 / (denCO2pump * epump);
Cpumpcaptot = Cpumpcapfix + Wpump * Cpumpcapvar;

cost = [Cpumpcaptot,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(464,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%O&M for Pressure Boosting Pump
cost_values(465,2) = operations(2);
cost_values(465,3) = operations(3);

pelec = 0.1036;
Cpumpelec = Wpump * 365*24 * pelec / maxCO2mult;
apumpOMvar	= 0.04;
CpumpOM = Cpumpcaptot * apumpOMvar;
CpumpOMtot = CpumpOM + Cpumpelec;

cost = [CpumpOMtot,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(465,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

% Distribution Pipeline Network

%Header Capital Cost
cost_values(466,2) = operations(2);
cost_values(466,3) = operations(2);

Nhead = ceil(Ninjwell/10);
Cheadcapvar	= 200000;
Cheadcaptot = Nhead * Cheadcapvar;

cost = [Cheadcaptot,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(466,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Header O&M
cost_values(467,2) = operations(2);
cost_values(467,3) = operations(3);

aheadOMvar = 0.04;
CheadOM = Cheadcaptot * aheadOMvar;

cost = [CheadOM,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(467,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Pipes to Inj wells Capital
cost_values(468,2) = operations(2);
cost_values(468,3) = operations(2);
mCO2actwell = mCO2actyr / Ninjwellactiv;

if mCO2actwell < 0.29
    Cpipwlcapvar = 360000;
elseif mCO2actwell < 0.676
    Cpipwlcapvar = 540000;
elseif mCO2actwell < 1.16
    Cpipwlcapvar = 720000;
else
    Cpipwlcapvar = 900000;
end

rCO2pl = sqrt(ACO2plz/ pi);
Lpipewell = apipewell * rCO2pl;
Cpipwlcaptot = Ninjwell * Lpipewell * Cpipwlcapvar;

cost = [Cpipwlcaptot,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(468,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Pipes to inj wells O&M
cost_values(469,2) = operations(2);
cost_values(469,3) = operations(3);

CpipwlOM = CpipfdOMvar * Lpipewell * Ninjwell;

cost = [CpipwlOM,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(469,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Building cost - Capital
cost_values(470,2) = operations(2);
cost_values(470,3) = operations(2);

cost = [200000,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(470,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Building O&M
cost_values(471,2) = operations(2);
cost_values(471,3) = operations(3);

cost = [0.04*200000,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(471,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


%Control System Capital
cost_values(472,2) = operations(2);
cost_values(472,3) = operations(2);

cost = [120000,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(472,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Control System O&M
cost_values(473,2) = operations(2);
cost_values(473,3) = operations(3);

cost = [0.04*120000,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(473,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Access Road to Building - Capital
cost_values(474,2) = operations(2);
cost_values(474,3) = operations(2);

RoadLen = 0.25; %mi
RoadCapvar = 1000000;%$/mi
Roadcap = RoadCapvar*RoadLen;

cost = [Roadcap,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(474,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Access Road to Building - O&M
cost_values(475,2) = operations(2);
cost_values(475,3) = operations(3);

aRoadOM = 0.04;
RoadOM = aRoadOM*Roadcap;

cost = [RoadOM,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(475,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Custody Transfer Gauge
cost_values(476,2) = 6;
cost_values(476,3) = 6;

cost = [250000,1,0,0,0,0,0,0,0,0,1,1,1,1+project_contingency];
cost_values(476,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%490 - ERR FEES for Carbon
cost_values(490,2) = operations(2);
cost_values(490,3) = operations(3);
cost_values(490,4) = 1;

ERR_Premium = 0.750; %ERR Premium($/tonne/yr):

Insurance = mCO2yr * ERR_Premium;
cost = [Insurance,1,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(490,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%% OTHER VARIABLE DEPENDENT COSTS
% Area dependent costs
cost_values(15:16,2) = site_selection_characterization(2);
cost_values(15:16,3) = site_selection_characterization(3);
cost_values(15:16,4) = 3;

cost = [10000,ACO2plUn,0,0,0,0,0,0,0,0,1,1,1+process_contingency,1+project_contingency];
cost_values(15,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [10000,ACO2plUn,0,0,0,0,0,0,0,0,1,1,1+process_contingency,1+project_contingency];
cost_values(16,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(83,2) = operations(2);
cost_values(83,3) = operations(3);
cost_values(83,4) = 5;

x = (PlumeUnArmult*PlumePrFAORmult*ACO2plz);
CAWells_total = round(x*dens_wells_CA);
CAwells_five = CAWells_total /(operations(1)/5);

cost = [31200,CAwells_five,11400,CAwells_five,13500,CAwells_five,0,0,0,0,1,1,1,1];
cost_values(83,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

%Fees for Injection

cost_values(84:86,2) = operations(2);
cost_values(84:86,3) = operations(3);
cost_values(84:86,4) = 1;

cost = [lifeTimeTonnes,0.25,1,0,0,0,0,0,0,0,1/operations(1),1,1,1];
cost_values(84,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [lifeTimeTonnes,0.07,1,0,0,0,0,0,0,0,1/operations(1),1,1,1];
cost_values(85,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost = [lifeTimeTonnes,0.01,1,0,0,0,0,0,0,0,1/operations(1),1,1,1];
cost_values(86,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(122,2) = operations(2);
cost_values(122,3) = operations(3);
cost_values(122,4) = 1;

cost = [lifeTimeTonnes/operations(1)/365/Ninjwell,0.05,0,0,0,0,0,0,0,0,1,1,1+process_contingency,1];
cost_values(122,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);



cost_values(131,2) = site_selection_characterization(2);
cost_values(131,3) = site_selection_characterization(3);
cost_values(131,4) = 1;

cost = [Ninjwell*20*10,4,25,4*20*Ninjwell,0,0,0,0,0,0,1,1,1+process_contingency,1];
cost_values(131,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


cost_values(132,2) = permitting_construction(2);
cost_values(132,3) = permitting_construction(3);
cost_values(132,4) = 1;

cost = [0,0,0,0,0,0,0,0,0,0,1,1,1,1];
cost_values(132,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);

cost_values(133,2) = operations(2);
cost_values(133,3) = operations(3);
cost_values(133,4) = 1;

cost = [Ninjwell*20*10,4,25,4*20*Ninjwell,0,0,0,0,0,0,1,1,1+process_contingency,1];
cost_values(133,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


cost_values(134,2) = pisc_closure(2);
cost_values(134,3) = pisc_closure(3);
cost_values(134,4) = 1;

cost = [Ninjwell*20*10,4,25,4*20*Ninjwell,0,0,0,0,0,0,1,1,1+process_contingency,1];
cost_values(134,1)=(cost(1)*cost(2)+cost(3)*cost(4)+cost(5)*cost(6)+cost(7)*cost(8)+cost(9)*cost(10))*cost(11)*cost(12)*cost(13)*cost(14);


end


%% Plume Well Schedule
function  [plume_area_other,monitoring_wells_reservoir,monitoring_wells_dual,monitoring_wells_above,monitoring_ground_water,monitoring_vadose] = plume_well_schedule(site_selection_characterization,permitting_construction,operations, pisc_closure,...
    PlumeUnArmult,ACO2plz,PlumePrFAORmult,Ang3DSeis,Ang2DSeis,frac_3D,frac_2D,ACO2plUn,Ltop,htot_ft,num_strat_test_wells,...
    strat_2_inj,Nwell_fin,well_spacing)

%{
Timelines:
Monitoring Wells	In Reservoir	7	36
	Above Seal	7	36
	Dual Completed	7	36
	Groundwater	7	36
	Vadose Zone	7	36

Injection Wells		6	6
Strat Test Wells		2	2
Water Production		7	36
Water Disposal		7	36
%}

strat=[site_selection_characterization(2),site_selection_characterization(2)];
inj = [permitting_construction(3),permitting_construction(3)];
monitoring = [operations(2),operations(3)];


Lbot = (Ltop + htot_ft)/5280; %feet converted to miles,depth to bottom of formation

CO2_Pressure_Front_AoR = (ACO2plUn*PlumePrFAORmult);

max_3D_seismic_A = pi*(sqrt(ACO2plUn/pi)+Lbot*tan(Ang3DSeis*pi/180))^2;%Maximum 3D Seismic Area

max_2D_seismic_L=2*(sqrt(ACO2plUn/pi)+Lbot*tan(Ang2DSeis*pi/180));%%Maximum 2D Seismic Length

start_A_3D_seis_monit =max_3D_seismic_A*frac_3D; %Starting area for 3D seismic monit.
start_L_2D_seis_monit =max_2D_seismic_L*frac_2D; %length for 2D seismic monit.


if strat_2_inj > Nwell_fin
    strat_2_inj_act = Nwell_fin;
else
    strat_2_inj_act = strat_2_inj;
end

strat_2_mon = sum(well_spacing(:,2))';

reg_strat_test_wells =num_strat_test_wells-(strat_2_inj+strat_2_mon);

reg_inj_well = Nwell_fin - strat_2_inj_act;%regular injection wells

num_water_production_wells = 0; %Water is considered off, Max. No. of Water Production Wells
num_water_disposal_wells = 0; %Water is off Max. No. of Water Disposal Wells

well_spacing(1,7) = ceil(ACO2plUn/well_spacing(1,1));
well_spacing(2,7) = ceil((CO2_Pressure_Front_AoR-ACO2plUn)/well_spacing(2,1));
well_spacing(3,7) = ceil(ACO2plUn/well_spacing(3,1));
well_spacing(4,7) = ceil((CO2_Pressure_Front_AoR-ACO2plUn)/well_spacing(4,1));
well_spacing(5,7) = ceil(ACO2plUn/well_spacing(5,1));
well_spacing(6,7) = ceil((CO2_Pressure_Front_AoR-ACO2plUn)/well_spacing(6,1));

well_spacing(:, 8) = max([well_spacing(:, 2), well_spacing(:, 3), well_spacing(:, 4)], [], 2);

well_spacing(:, 9) = min([max([well_spacing(:, 5), zeros(size(well_spacing, 1), 1) ], [], 2),max([well_spacing(:, 7), well_spacing(:, 8)], [], 2)],[],2);

well_spacing(:,10) = max([well_spacing(:, 2), zeros(size(well_spacing, 1), 1),well_spacing(:, 9)], [], 2);

%% Plume Flows
%CO2 Plume Area & Other Area Calcs
plume_area_other = zeros(24, pisc_closure(3));  % Preallocate the plume_area_other array

i = 1:pisc_closure(3);

%CO2 Plume
plume_area_other(1,i) = ((i >= operations(2)) .* min((ACO2plz / operations(1)) * (i - operations(2) + 1), ACO2plz)) + ((i < operations(2)) .* 0);
plume_area_other(2,i) = (i < operations(2)) .* (ACO2plz * PlumeUnArmult) + (i >= operations(2)) .* (plume_area_other(1,i) * PlumeUnArmult);
plume_area_other(3,i) = plume_area_other(2,i) .* PlumePrFAORmult;
plume_area_other(4,i) = (i < operations(2)) .* plume_area_other(2,i) + (i >= operations(2)) .* (i <= operations(3)) .* (ACO2plUn * (frac_3D + (1 - frac_3D) .* (i - operations(2)) / (operations(3) - operations(2)))) + (i > operations(3)) .*ACO2plUn;
plume_area_other(5,i) = pi*(sqrt(plume_area_other(4,i)/pi)+Lbot*tan(Ang3DSeis*pi/180)).^2;
plume_area_other(6,i) = (i < operations(2)) .* plume_area_other(2,i) + (i >= operations(2)) .* (i <= operations(3)) .* (ACO2plUn * (frac_2D + (1 - frac_2D) .* (i - operations(2)) / (operations(3) - operations(2)))) + (i > operations(3)) .*ACO2plUn;
plume_area_other(7,i) = 2*(sqrt(plume_area_other(6,i)/pi)+Lbot*tan(Ang2DSeis*pi/180));

%All Types of Strat Test Wells
plume_area_other(9,i) = (i >= strat(1)).*num_strat_test_wells + (i < strat(1)) .* 0;
plume_area_other(8, i) = plume_area_other(9, i) .* (i == 1) + (plume_area_other(9, i) - plume_area_other(9, max(i - 1, 1))) .* (i ~= 1);

%Regular Strat Test Wells
plume_area_other(11,i) = (i >= strat(1)).*reg_strat_test_wells + (i < strat(1)) .* 0;
plume_area_other(10, i) = plume_area_other(11, i) .* (i == 1) + (plume_area_other(11, i) - plume_area_other(11, max(i - 1, 1))) .* (i ~= 1);

%Strat Wells to Mon. Wells.
plume_area_other(13,i) = (i >= strat(1)).*strat_2_mon + (i < strat(1)) .* 0;
plume_area_other(12, i) = plume_area_other(13, i) .* (i == 1) + (plume_area_other(13, i) - plume_area_other(13, max(i - 1, 1))) .* (i ~= 1);

%Injection Wells
%Strat wells converted to inj wells
plume_area_other(14,i) = (i == strat(1)).*strat_2_inj_act + (i ~= strat(1)) .* 0;
plume_area_other(15,i) = (i == inj(1)).*strat_2_inj_act + (i ~= inj(1)) .* 0;
plume_area_other(16,i) = (i >= inj(1)).*strat_2_inj_act + (i < inj(1)) .* 0;

%Regular Injection wells
plume_area_other(18,i) = (i >= inj(1)).*reg_inj_well + (i < inj(1)) .* 0;
plume_area_other(17, i) = plume_area_other(18, i) .* (i == 1) + (plume_area_other(18, i) - plume_area_other(18, max(i - 1, 1))) .* (i ~= 1);

plume_area_other(19,i) = plume_area_other(17, i) + plume_area_other(15, i);
plume_area_other(20,i) = plume_area_other(18, i) + plume_area_other(16, i);

plume_area_other(22,i) = (i >= monitoring(1)).*num_water_production_wells+ (i < monitoring(1)) .* 0;
plume_area_other(21, i) = plume_area_other(22, i) .* (i == 1) + (plume_area_other(22, i) - plume_area_other(22, max(i - 1, 1))) .* (i ~= 1);

plume_area_other(24,i) = (i >= monitoring(1)).*num_water_disposal_wells + (i < monitoring(1)) .* 0;
plume_area_other(23, i) = plume_area_other(24, i) .* (i == 1) + (plume_area_other(24, i) - plume_area_other(24, max(i - 1, 1))) .* (i ~= 1);

%Monitoring Wells Schedules
monitoring_wells_reservoir = zeros(13, pisc_closure(3));

monitoring_wells_reservoir(3,i) = (i >= monitoring(1)).*ceil(plume_area_other(2,i).*(well_spacing(1,10)/ACO2plUn)) + (i < monitoring(1)) .* 0;
monitoring_wells_reservoir(4,i) = (i == monitoring(1)).*well_spacing(1,2) + (i ~= monitoring(1)).*0;
monitoring_wells_reservoir(6,i) = (well_spacing(1,10) == 0) * 0 + (well_spacing(1,10) ~= 0) .* (i < monitoring(1)) .* 0 + (well_spacing(1,10) ~= 0) .* (i == monitoring(1)) .* max(well_spacing(1,6), monitoring_wells_reservoir(3,i)) + (well_spacing(1,10) ~= 0) * (i > monitoring(1)) .* max(monitoring_wells_reservoir(3,i), monitoring_wells_reservoir(6,max(i-1,1)));
monitoring_wells_reservoir(5,i) = (i == 1) .* monitoring_wells_reservoir(6,i) + (i ~= 1) .* monitoring_wells_reservoir(6,i) - monitoring_wells_reservoir(6,max(i-1,1)) - monitoring_wells_reservoir(4,i);

monitoring_wells_reservoir(1,i) = ((monitoring_wells_reservoir(6,i) > 0) .* plume_area_other(2,i)./max(1,monitoring_wells_reservoir(6,i))) + (monitoring_wells_reservoir(6,i) <= 0) * 0;

monitoring_wells_reservoir(7,i) = (i >= monitoring(1)).*ceil((plume_area_other(3,i)-plume_area_other(2,i)).*(well_spacing(2,10)/(CO2_Pressure_Front_AoR-ACO2plUn))) + (i < monitoring(1)) .* 0;
monitoring_wells_reservoir(8,i) = (i == monitoring(1)).*well_spacing(2,2) + (i ~= monitoring(1)).*0;
monitoring_wells_reservoir(10,i) = (well_spacing(2,10) == 0) * 0 + (well_spacing(2,10) ~= 0) .* (i < monitoring(1)) .* 0 + (well_spacing(2,10) ~= 0) .* (i == monitoring(1)) .* max(well_spacing(2,6), monitoring_wells_reservoir(7,i)) + (well_spacing(2,10) ~= 0) * (i > monitoring(1)) .* max(monitoring_wells_reservoir(7,i), monitoring_wells_reservoir(10,max(i-1,1)));
monitoring_wells_reservoir(9,i) = (i == 1) .* monitoring_wells_reservoir(10,i) + (i ~= 1) .* monitoring_wells_reservoir(10,i) - monitoring_wells_reservoir(10,max(i-1,1)) - monitoring_wells_reservoir(8,i);

monitoring_wells_reservoir(2,i) = ((monitoring_wells_reservoir(10,i) > 0) .* (plume_area_other(3,i)-plume_area_other(2,i))./max(1,monitoring_wells_reservoir(10,i))) + (monitoring_wells_reservoir(10,i) <= 0) * 0;

monitoring_wells_reservoir(11,i) = monitoring_wells_reservoir(4,i)+monitoring_wells_reservoir(8,i);
monitoring_wells_reservoir(12,i) = monitoring_wells_reservoir(5,i)+monitoring_wells_reservoir(9,i);
monitoring_wells_reservoir(13,i) = monitoring_wells_reservoir(6,i)+monitoring_wells_reservoir(10,i);

%Monitoring Wells above Seal
monitoring_wells_above = zeros(13, pisc_closure(3));

monitoring_wells_above(3,i) = (i >= monitoring(1)).*ceil(plume_area_other(2,i).*(well_spacing(3,10)/ACO2plUn)) + (i < monitoring(1)) .* 0;
monitoring_wells_above(4,i) = (i == monitoring(1)).*well_spacing(3,2) + (i ~= monitoring(1)).*0;
monitoring_wells_above(6,i) = (well_spacing(3,10) == 0) * 0 + (well_spacing(3,10) ~= 0) .* (i < monitoring(1)) .* 0 + (well_spacing(3,10) ~= 0) .* (i == monitoring(1)) .* max(well_spacing(3,6), monitoring_wells_above(3,i)) + (well_spacing(3,10) ~= 0) * (i > monitoring(1)) .* max(monitoring_wells_above(3,i), monitoring_wells_above(6,max(i-1,1)));
monitoring_wells_above(5,i) = (i == 1) .* monitoring_wells_above(6,i) + (i ~= 1) .* monitoring_wells_above(6,i) - monitoring_wells_above(6,max(i-1,1)) - monitoring_wells_above(4,i);

monitoring_wells_above(1,i) = ((monitoring_wells_above(6,i) > 0) .* plume_area_other(2,i)./max(1,monitoring_wells_above(6,i))) + (monitoring_wells_above(6,i) <= 0) * 0;

monitoring_wells_above(7,i) = (i >= monitoring(1)).*ceil((plume_area_other(3,i)-plume_area_other(2,i)).*(well_spacing(4,10)/(CO2_Pressure_Front_AoR-ACO2plUn))) + (i < monitoring(1)) .* 0;
monitoring_wells_above(8,i) = (i == monitoring(1)).*well_spacing(4,2) + (i ~= monitoring(1)).*0;
monitoring_wells_above(10,i) = (well_spacing(4,10) == 0) * 0 + (well_spacing(4,10) ~= 0) .* (i < monitoring(1)) .* 0 + (well_spacing(4,10) ~= 0) .* (i == monitoring(1)) .* max(well_spacing(4,6), monitoring_wells_above(7,i)) + (well_spacing(4,10) ~= 0) * (i > monitoring(1)) .* max(monitoring_wells_above(7,i), monitoring_wells_above(10,max(i-1,1)));
monitoring_wells_above(9,i) = (i == 1) .* monitoring_wells_above(10,i) + (i ~= 1) .* monitoring_wells_above(10,i) - monitoring_wells_above(10,max(i-1,1)) - monitoring_wells_above(8,i);

monitoring_wells_above(2,i) = ((monitoring_wells_above(10,i) > 0) .* (plume_area_other(3,i)-plume_area_other(2,i))./max(1,monitoring_wells_above(10,i))) + (monitoring_wells_above(10,i) <= 0) * 0;

monitoring_wells_above(11,i) = monitoring_wells_above(4,i)+monitoring_wells_above(8,i);
monitoring_wells_above(12,i) = monitoring_wells_above(5,i)+monitoring_wells_above(9,i);
monitoring_wells_above(13,i) = monitoring_wells_above(6,i)+monitoring_wells_above(10,i);

%Monitoring Wells Dual
monitoring_wells_dual = zeros(13, pisc_closure(3));

monitoring_wells_dual(3,i) = (i >= monitoring(1)).*ceil(plume_area_other(2,i).*(well_spacing(5,10)/ACO2plUn)) + (i < monitoring(1)) .* 0;
monitoring_wells_dual(4,i) = (i == monitoring(1)).*well_spacing(5,2) + (i ~= monitoring(1)).*0;
monitoring_wells_dual(6,i) = (well_spacing(5,10) == 0) * 0 + (well_spacing(5,10) ~= 0) .* (i < monitoring(1)) .* 0 + (well_spacing(5,10) ~= 0) .* (i == monitoring(1)) .* max(well_spacing(5,6), monitoring_wells_dual(3,i)) + (well_spacing(5,10) ~= 0) * (i > monitoring(1)) .* max(monitoring_wells_dual(3,i), monitoring_wells_dual(6,max(i-1,1)));
monitoring_wells_dual(5,i) = (i == 1) .* monitoring_wells_dual(6,i) + (i ~= 1) .* monitoring_wells_dual(6,i) - monitoring_wells_dual(6,max(i-1,1)) - monitoring_wells_dual(4,i);

monitoring_wells_dual(1,i) = ((monitoring_wells_dual(6,i) > 0) .* plume_area_other(2,i)./max(1,monitoring_wells_dual(6,i))) + (monitoring_wells_dual(6,i) <= 0) * 0;

monitoring_wells_dual(7,i) = (i >= monitoring(1)).*ceil((plume_area_other(3,i)-plume_area_other(2,i)).*(well_spacing(6,10)/(CO2_Pressure_Front_AoR-ACO2plUn))) + (i < monitoring(1)) .* 0;
monitoring_wells_dual(8,i) = (i == monitoring(1)).*well_spacing(6,2) + (i ~= monitoring(1)).*0;
monitoring_wells_dual(10,i) = (well_spacing(6,10) == 0) * 0 + (well_spacing(6,10) ~= 0) .* (i < monitoring(1)) .* 0 + (well_spacing(6,10) ~= 0) .* (i == monitoring(1)) .* max(well_spacing(6,6), monitoring_wells_dual(7,i)) + (well_spacing(6,10) ~= 0) * (i > monitoring(1)) .* max(monitoring_wells_dual(7,i), monitoring_wells_dual(10,max(i-1,1)));
monitoring_wells_dual(9,i) = (i == 1) .* monitoring_wells_dual(10,i) + (i ~= 1) .* monitoring_wells_dual(10,i) - monitoring_wells_dual(10,max(i-1,1)) - monitoring_wells_dual(8,i);

monitoring_wells_dual(2,i) = ((monitoring_wells_dual(10,i) > 0) .* (plume_area_other(3,i)-plume_area_other(2,i))./max(1,monitoring_wells_dual(10,i))) + (monitoring_wells_dual(10,i) <= 0) * 0;

monitoring_wells_dual(11,i) = monitoring_wells_dual(4,i)+monitoring_wells_dual(8,i);
monitoring_wells_dual(12,i) = monitoring_wells_dual(5,i)+monitoring_wells_dual(9,i);
monitoring_wells_dual(13,i) = monitoring_wells_dual(6,i)+monitoring_wells_dual(10,i);

%Monitoring Wells Groundwater
monitoring_ground_water = zeros(2, pisc_closure(3));
monitoring_ground_water(2,i) = (i >= monitoring(1)) .* well_spacing(7,1) .* plume_area_other(20,i) + (i < monitoring(1)) .* 0;
monitoring_ground_water(1,i) = (i == 1) .* monitoring_ground_water(2,i) + (i ~= 1) .* monitoring_ground_water(2,i) - monitoring_ground_water(2,max(i-1,1));

%Monitoring Wells Vadose Zone
monitoring_vadose= zeros(2, pisc_closure(3));
monitoring_vadose(2,i) = (i >= monitoring(1)) .* well_spacing(8,1) .* plume_area_other(20,i) + (i < monitoring(1)) .* 0;
monitoring_vadose(1,i) = (i == 1) .* monitoring_vadose(2,i) + (i ~= 1) .* monitoring_vadose(2,i) - monitoring_vadose(2,max(i-1,1));

end



