function [transport_costs,dpt]=transportModule(tonnes,lengthMi,capFact,escalation_rate,discount_rates)
    format("bank")
    
    load("Transport_Data.mat")
    op_years = 12;

    numPumps = calcOptPumps(tonnes,lengthMi,capFact);
    
    
    
    minInDia=dia_in_min(tonnes,lengthMi,capFact,numPumps);
    
    if isstring(minInDia)
        transport_costs = inf;
        dpt = [inf;inf;inf];
        opex = inf;
        capex = inf;
        capexNom = inf;
        opexNom = inf;
    else
    nomPipeDiam=Pipe_Size(minInDia);

    %pipeline costs from Parker Eqn
    capex=parkerEq(lengthMi, nomPipeDiam);
    qm = (tonnes/capFact)/1000000; %convert from tonnes to megatonnes per year
    Dens_CO2 = denDuanCO2(0, AveragePressure, temp);
    pumpPower =(1000000000*1000000*0.001/(365*24*3600))*qm*dPressure/(pumpEfficency*Dens_CO2); %in KW
    %additional pipeline costs in 2011 $
    capex.surgeTank_Cost = 1244744;
    capex.pipelineControlSystem_Cost = 111907;
    capex.total = capex.total + capex.surgeTank_Cost + capex.pipelineControlSystem_Cost;
    capex.total_pipeline = capex.total;
    capex.pump_cost = (fixed_pump_cost + variable_pump_cost * pumpPower) * numPumps * convert2005to2011;
    capex.total = capex.total + capex.pump_cost;

    %opex pipeline O&M
    opex.pipeline = capex.total_pipeline * O_M_percent_pipeline;
    opex.equipment = (capex.pump_cost + capex.surgeTank_Cost + capex.pipelineControlSystem_Cost) * O_M_percent_equipment;
    electricityPerYear = capFact * 8760 * numPumps * pumpPower/1000; %MWhr/yr
    opex.electricity = electricityPerYear * electricityCost;
    opex.total = opex.pipeline + opex.equipment + opex.electricity;

    
    capexNom = zeros(1,op_years+duration_of_construction);
    opexNom = zeros(1,op_years+duration_of_construction);

    esc_sched = zeros(1,op_years+duration_of_construction);
    disc_sched_3 = zeros(1,op_years+duration_of_construction);
    disc_sched_2 = zeros(1,op_years+duration_of_construction);
    
    diff = start_year-base_year;
    for i=1:(op_years+duration_of_construction)
        esc_sched(i) = (1+escalation_rate)^(diff+i-1);
        disc_sched_3(i) = 1/(1+discount_rates(1,1))^(i-1);
        disc_sched_2(i) = 1/(1+discount_rates(1,2))^(i-1);
        if i<= 3
           capexNom(i) = capex.material_costs * brMat(i) + capex.labor_costs*brLabor(i) + capex.rightofway_costs*brROW(i) + ...
            capex.miscellanous_costs*brMisc(i) + capex.surgeTank_Cost*brSurge(i) + capex.pipelineControlSystem_Cost*brControl(i) + capex.pump_cost*brPump(i);
        else
            opexNom(i) = opex.total;
        end
    end
    capexNom = capexNom.*(1+contingency) .* esc_sched;
    opexNom = opexNom .* esc_sched;

    PV_capex_3 = capexNom.*disc_sched_3;
    PV_opex_3 = opexNom.*disc_sched_3;
    
    PV_capex_2 = capexNom.*disc_sched_2;
    PV_opex_2 = opexNom.*disc_sched_2;
    
    transport_costs=[capexNom;opexNom;PV_capex_3;PV_opex_3;PV_capex_2;PV_opex_2];

    dpt_nom = (sum(capexNom) + sum(opexNom))/(op_years*tonnes);
    dpt_PV_3 = (sum(PV_capex_3) + sum(PV_opex_3))/(op_years*tonnes);
    dpt_PV_2 = (sum(PV_capex_2) + sum(PV_opex_2))/(op_years*tonnes);
    
    dpt = [dpt_nom;dpt_PV_3;dpt_PV_2];
    end

end

function capex=parkerEq(lengthMi, diam)

if diam <= 12
    fnc = 1;
elseif diam <= 16
    fnc = 1.12;
elseif diam <= 20
    fnc = 1.18;
else
    fnc = 1.25;
end

%Natural gas Pipeline cost
% $2000
Cmat1 = (35000 + lengthMi * (330.5 * diam ^ 2 + 687 * diam + 26960))/lengthMi;
Clab1 = (185000 + lengthMi * (343 * diam ^ 2 + 2074 * diam + 170013))/lengthMi;
Crow1 = (40000 + lengthMi * (577 * diam + 29788))/lengthMi;
Cmisc1 = (95000 + lengthMi * (8417 * diam + 7324))/lengthMi;

%converting to $2011
Cmat = Cmat1 * (525 / 261.0); %'Handy-Whitman index for gas transmission pipelines
Clab = Clab1 * (525 / 261.0); %'Handy-Whitman index for gas transmission pipelines
Crow = Crow1 * (113.8 / 88.7); %'Gross domestic product chain-type price index
Cmisc = Cmisc1 * (190.9 / 122.3); %'Producer price index

%CO2 pipeline
capex.material_costs_pmi = Cmat *fnc;
capex.labor_costs_pmi = Clab *fnc;
capex.rightofway_costs_pmi = Crow;
capex.miscellanous_costs_pmi = Cmisc;
%cost per mile
capex.material_costs = Cmat *fnc * lengthMi;
capex.labor_costs = Clab *fnc * lengthMi;
capex.rightofway_costs = Crow * lengthMi;
capex.miscellanous_costs = Cmisc * lengthMi;
COtot_mi = Cmat*fnc + Clab*fnc + Crow + Cmisc;
capex.total = COtot_mi * lengthMi;
end

function numPumps = calcOptPumps(tonnes,lengthMi,capFact)
   
    Dia_nom = 48; %finds the maximum length of a pipe segment when the nominal pipe size is 48 inches, the largest pipe size currently in the model 
    Dia_in = Dia_in_nom(Dia_nom);
    l_max = L_seg_max(Dia_in,tonnes,lengthMi,capFact);
    l_max_mi = l_max * 0.000621371;

    if l_max_mi > lengthMi
        Npump_min = 0;
        minInDia=dia_in_min(tonnes,lengthMi,capFact,Npump_min);
        if isstring(minInDia)
        else
            Dia_nom=Pipe_Size(minInDia);
        end
    else
        Dia_nom = 48;
        Nseg_1 = lengthMi / l_max_mi;
        Nseg = ceil(Nseg_1);
        Npump = Nseg - 1;
        Npump_min = Npump;
    end
  

%     Npump_min = 0;
%     pr_CO2_min = 1e+99;
%     for idst = 0:20
%         [~,~,~,~,pr_pipe] = dollarPerTonne(idst,tonnes,lengthMi,capFact);
%         if pr_pipe < pr_CO2_min 
%            Npump_min = idst;
%            pr_CO2_min = pr_pipe;
%         end
% 
%     end
% 

    %{
    pr_CO2_min = 1e+99;

    for idst=1:100
        Dia_in = Dia_in_nom(Dia_nom);
        l_max = L_seg_max(Dia_in,tonnes,lengthMi,capFact);
        l_max_mi = l_max * 0.000621371;
        Nseg_1 = lengthMi / l_max_mi;
        if Nseg_1 > 200 * (Npump_min +1)
            break;
        else
            Nseg = ceil(Nseg_1);
            Npump = Nseg - 1;
        end
        
        
       [~,~,~,~,pr_pipe] = dollarPerTonne(Npump,tonnes,lengthMi,capFact);
        if  pr_pipe < pr_CO2_min
            pr_CO2_min = pr_pipe;
            Npump_min = Npump;
        end
        Dia_nom = Dia_next_low(Dia_nom);
        
        if Dia_nom < 0.0001
            break;
        end
        
    end 
    %}
    numPumps = Npump_min;
    
end

%function returns the inner diameter in inches of a standard pipe with a
%nominal diameter 

function PS = Dia_in_nom(dia)
Dia_nom = dia;

if Dia_nom <= 4
    PS = 4;
elseif Dia_nom <= 6
    PS = 6;
elseif Dia_nom <= 8
    PS = 8;
elseif Dia_nom <= 10
    PS = 10;
elseif Dia_nom <= 12
    PS = 12;
elseif Dia_nom <= 16
    PS = 15.162;
elseif Dia_nom <= 20
    PS = 18.952;
elseif Dia_nom <= 24
    PS = 22.742;
elseif Dia_nom <= 30
    PS = 28.428;
elseif Dia_nom <= 36
    PS =34.114;
elseif Dia_nom <= 42
    PS =39.8;
elseif Dia_nom <= 48
    PS = 45.5;
else
    PS = 99.9;
end

end 


%converts minimum pipe diameter to nominal or standard diameter of
%commercially availble pipe
%returns nominal pipeline diameter
function PS=Pipe_Size(dia)

%converting meters to inches
Dia_in = dia * 39.3701;

if Dia_in <= 4 
    PS = 4;
elseif Dia_in <= 6 
    PS = 6;
elseif Dia_in <= 8
    PS = 8;
elseif Dia_in <= 10
    PS = 10;
elseif Dia_in <= 12
    PS = 12;
elseif Dia_in <= 15.162 
    PS = 16;
elseif Dia_in <= 18.952 
    PS = 20;
elseif Dia_in <= 22.742
    PS = 24;
elseif Dia_in <= 28.428 
    PS = 30;
elseif Dia_in <= 34.114
    PS = 36;
elseif Dia_in <= 39.8
    PS = 42;
elseif Dia_in <= 45.5
    PS = 48;
else
    PS = 2000;
end

end

%finds the minimum pipe inner diameter that can sustain the flow rate given
%by q_max over the pressure change given by p_in, p_out and accounts for
%elevation change across pipeline, outputs minInDia in meters
function minInDia=dia_in_min(tonnes,lengthMi,capFact,numPumps)

g = 9.80665; %acceleration due to gravity

lengthM = lengthMi/ 0.621371 * 1000;
elevation = 0; % ASSUMING 0 ELEVATION CHANGE
etaIn = 0.0018; %Pipe inside surface roughness in inches
etaMe = etaIn * 0.0254; %converting to meters

p_in = 2200 * 6894.76;%converting from psig to pa
p_out = 1200 * 6894.76;%converting from psig to pa

q_max = (tonnes/capFact)*1000/(365*24*60*60); %maximum CO2 flow rate converted to kg/s
temp = (53 - 32) * 5.0/9.0 + 273.15; %converting temp in to Kelvin

Nseg = numPumps + 1; %calculates the number of pipe segments 

L_seg = lengthM/Nseg; %calculates the length of a segment
h_dif_seg = elevation/Nseg; %elevation increase or decrease along the segment


%Method to find Min Inner Diamerter for incompressible fluid

p_liq_av = 0.5 * (p_out + p_in);
Dens_CO2 = denDuanCO2(5E-05, p_liq_av * 1E-06, temp);
Visc_CO2 = visFWVCO2(Dens_CO2, temp) * 1E-06;

Dia_g = 0.5;%variables used for iteration, in meters
Dia_old = 0.9 * Dia_g;

%initalize variables used to check for denominator that leades to complex number for a result
denom_chk = 0;

%initialize counter
count =0;

while ~(count > 1000 || abs(Dia_old/Dia_g - 1) < 10^-6 || denom_chk  == 1)
    %Calculate Fanning friction factor

    %calculate Reynolds number for SC/Liquid CO2
    Re = 4*q_max/(pi*Visc_CO2*Dia_g);
    %Colebrook-White fanning friction factor
    ff = FFF_Cole(Re,Dia_g,etaMe);

    %In denominator, p_in - p_out i spressure difference from outlet to
    %inlet. h_dif_seg is elevation increase or decrease from inlet to
    %outlet (h_out - h_in_ of the pipe segment, the formula requires h_in -
    %h_out
    num = 32 * ff * q_max^2 * L_seg;
    denom = pi^2 * Dens_CO2 * ((p_in - p_out) - g * Dens_CO2 * h_dif_seg);
    if denom <= 0
        % Since denom is not negative then the quantity
        % num / denom will be negative and when raised to the 1/5 power
        % will give a result that is an imaginary number which is not
        % physically realistic.
        % A result that is an imaginary number implies that the
        % specified pressure drop can never overcome the
        % frictional losses in the pipe and also the potential
        % energy due to an elevation increase if the pipe segment has an
        % elevation increase.
        % A shorter pipe segment is needed.
        % The algorithm returns a diameter of 99.9 inches (converted to meters)
        % to indicate this situation.
        % Also denom_chk is set to 1 to exit the loop.
        Dia_liq = 99.9/39.3701;
        denom_chk = 1;
    else
        %denom is negative which it must be since num is negative
        %otherwise the quantity num / denom yields a negative number
        %which when raised to the 0.2 power yields an imaginary number
        %which is not physically realistic
        Dia_liq = (num/denom)^0.2;
    end

    %update pipeline diameter
    Dia_old = Dia_g;
    Dia_g = Dia_liq;

    count = count+1;
end

if count > 1000 
    %no convergence
    Dia_liq = "No convergence";
end

%return value
minInDia = Dia_liq;

end

%calculates mass density from the molar volume resulting from the equation
%of state for CO2
function Dens_CO2=denDuanCO2(volinit,pres,tmp)

MW = 44.0095;
Dens_CO2 = MW * 0.001 / volDuanCO2(volinit,pres,tmp);
end


function result=volDuanCO2(volinit,pres,tmp)

Rgas_MPa = 0.000008314; %m3-MPa/K-mol
tol = 0.00001;

tmpcrit = 304.1282;%critical temperature
prescrit = 7.3773;%critical pressure (MPa)
dencrit = 467.6;
MW = 44.0095;
volcrit = MW * 0.001 / dencrit;

if volinit <=  0
    volla = Rgas_MPa * tmp / pres;

    if tmp > tmpcrit
        if pres > prescrit
            volla = volla * 1.05;
            if tmp > 500 && pres > 30
                volla = volla * 5;
            end
        else
            volla = volla * 0.95;
        end
    else

        %constants
        a1 = -7.0602087;
        a2 = 1.9391218;
        a3 = -1.6463597;
        a4 = -3.2995634;
        n1 = 1;
        n2 = 1.5;
        n3 = 2;
        n4 = 4;

        %Calculate vapor pressure
        vpSWCO2 = prescrit * exp((tmpcrit / tmp) * (a1 * (1 - tmp / tmpcrit) ^ n1 + a2 * (1 - tmp / tmpcrit) ^ n2 + a3 * (1 - tmp / tmpcrit) ^ n3 + a4 * (1- tmp / tmpcrit) ^ n4));
        pvap = vpSWCO2;
        if pres > pvap
            % Constants
            a1 = 1.9245108;
            a2 = -0.62385555;
            a3 = -0.32731127;
            a4 = 0.39245142;
            n1 = 0.34;
            n2 = 0.5;
            n3 = 10.0 / 6.0;
            n4 = 11.0 / 6.0;

            %Calculate vapor pressure
            denslSWCO2 = dencrit * exp(a1 * (1 - tmp / tmpcrit) ^ n1 + a2 * (1 - tmp / tmpcrit) ^ n2 + a3 * (1 - tmp / tmpcrit) ^ n3 + a4 * (1 - tmp / tmpcrit) ^ n4);
            vollsat = MW * 0.001 / denslSWCO2;
            if volla > 0.95 * vollsat
                volla = 0.95 * vollsat;
            end
        else
            % Constants
            a1 = -1.7074879;
            a2 = -0.8227467;
            a3 = -4.6008549;
            a4 = -10.111178;
            a5 = -29.742252;
            n1 = 0.34;
            n2 = 0.5;
            n3 = 1;
            n4 = 7.0 / 3.0;
            n5 = 14.0 / 3.0;
            %Calculate vapor pressure
            densvSWCO2 = dencrit * exp(a1 * (1 - tmp / tmpcrit) ^ n1 + a2 * (1 - tmp / tmpcrit) ^ n2 + a3 * (1 - tmp / tmpcrit) ^ n3 + a4 * (1 - tmp / tmpcrit) ^ n4 + a5 * (1 - tmp / tmpcrit) ^ n5);
            volvsat = MW * 0.001 /densvSWCO2;
            if volla < 1.05 * volvsat
                volla = 1.05*volvsat;
            end
        end
    end
else
    volla = volcrit;
end

for i=1:100
    if i == 1
    else
        if vol2 > 2
            break
        else
            volla = volla * 1.05;
        end
    end
    voll = volla;
    for j=1: 1000
        vol2 = voll - fDuanCO2(voll,pres,tmp)/d1fDuanCO2(voll,pres,tmp);
        reldif = abs(voll - vol2)/abs(vol2);
        if(reldif <= tol)
            break
        else
            voll = vol2;
        end
    end
end
result = vol2;

end

%calculates the nonlinear function resulting from the equation of state for
%CO2 developed by Duan et al

function result=fDuanCO2(vol,pres,tmp)
    
Rgas_Lbar = 0.08314467; % universal gas constant (L-bar/K-mol)

% Properties of CO2
tmpcrit = 304.1282; %critical temperature (K)
prescrit = 7.3773; %critical pressure (MPa)
volcrit = Rgas_Lbar * tmpcrit / (prescrit * 10); %not really critical volume (L/mol)

%Constants for Duan et al. CO2 equation of state
a1 = 0.0899288497;
a2 = -0.494783127;
a3 = 0.0477922245;
a4 = 0.0103808883;
a5 = -0.0282516861;
a6 = 0.0949887563;
a7 = 0.00052060088;
a8 = -0.000293540971;
a9 = -0.00177265112;
a10 = -2.51101973E-05;
a11 = 8.93353441E-05;
a12 = 7.88998563E-05;
a13 = -0.0166727022;
a14 = 1.398;
a15 = 0.0296;

%Calculate Tr, Pr, and Vr
Tr = tmp / tmpcrit;
Pr = pres / prescrit;
Vr = vol * 1000/ volcrit;

%Calculate intermediate values
b1 = a1 + (a2 / Tr ^ 2) + (a3 / Tr ^ 3);
b2 = a4 + (a5 / Tr ^ 2) + (a6 / Tr ^ 3);
b3 = a7 + (a8 / Tr ^ 2) + (a9 / Tr ^ 3);
b4 = a10 + (a11 / Tr ^ 2) + (a12 / Tr ^ 3);
b5 = a13 / Tr ^ 3;

%Calculate function
result = Pr * Vr ^ 6 / Tr - Vr ^ 5 - b1 * Vr ^ 4 - b2 * Vr ^ 3 - b3 * Vr - b4 - b5 * Vr * (a14 * Vr ^ 2 + a15) * exp(-a15 / Vr ^ 2);

end

%calculates the first derivative with respect to vol of the nonlinear
%function
function result = d1fDuanCO2(vol, pres, tmp)
    
Rgas_Lbar = 0.08314467;  %universal gas constant (L-bar/K-mol)

%Properties of CO2
tmpcrit = 304.1282; %critical temperature (K)
prescrit = 7.3773;  %critical pressure (MPa)
volcrit = Rgas_Lbar * tmpcrit / (prescrit * 10); % not really critical volume (L/mol)

%Constants for Duan et al. CO2 equation of state
a1 = 0.0899288497;
a2 = -0.494783127;
a3 = 0.0477922245;
a4 = 0.0103808883;
a5 = -0.0282516861;
a6 = 0.0949887563;
a7 = 0.00052060088;
a8 = -0.000293540971;
a9 = -0.00177265112;
%a10 = -2.51101973E-05;
%a11 = 8.93353441E-05;
%a12 = 7.88998563E-05;
a13 = -0.0166727022;
a14 = 1.398;
a15 = 0.0296;

%Calculate Tr, Pr, and Vr
Tr = tmp / tmpcrit;
Pr = pres / prescrit;
Vr = vol * 1000 / volcrit;

%Calculate intermediate values
b1 = a1 + (a2 / Tr ^ 2) + (a3 / Tr ^ 3);
b2 = a4 + (a5 / Tr ^ 2) + (a6 / Tr ^ 3);
b3 = a7 + (a8 / Tr ^ 2) + (a9 / Tr ^ 3);
%b4 = a10 + (a11 / Tr ^ 2) + (a12 / Tr ^ 3);
b5 = a13 / Tr ^ 3;

%Calculate d f / d Vr
dfdVr = 6 * Pr * Vr ^ 5 / Tr - 5 * Vr ^ 4 - 4 * b1 * Vr ^ 3 - 3 * b2 * Vr ^ 2 - b3 - (3 * b5 * a14 * Vr ^ 2 + b5 * a15) * exp(-a15 / Vr ^ 2) - (b5 * a14 * Vr ^ 3 + b5 * a15 * Vr) * (2 * a15 / Vr ^ 3) * exp(-a15 / Vr ^ 2);

%Need to return d f / d vol
%Recall that Vr = vol * 1000 / volcrit
%Thus, d f / d Vr = d f / d (vol * 1000 / volcrit)
%Since 1000 and volcrit are constants, we get
% d f / d Vr = (volcrit / 1000) d f / d vol Or
%d f / d vol = (1000 / volcrit) d f / d Vr
result = (1000 / volcrit) * dfdVr;

end

function Visc_CO2=visFWVCO2(den, tmp)

%Calculation of zero density viscosity term
tmpcs = 251.196;  %Kelvin
az(1) = 0.235156;
az(2) = -0.491266;
az(3) = 0.05211155;
az(4) = 0.05347906;
az(5) = -0.01537102;

CSR = 0;
for j=1:5
    CSR = CSR + az(j) * (log(tmp / tmpcs)) ^ (j - 1);
end

visz = 1.00697 * sqrt(tmp) / exp(CSR);

%Calculation of excess viscosity term that corrects for non-zero density
d11 = 0.004071119;
d21 = 7.198037E-05;
d64 = 2.411697E-17;
d81 = 2.971072E-23;
d82 = -1.627888E-23;

visexc = d11 * den + d21 * den ^ 2 + d64 * den ^ 6 * (tmp / tmpcs) ^ -3 + d81 * den ^ 8 + d82 * den ^ 8 * (tmp / tmpcs) ^ -1;

%Calculate viscosity
Visc_CO2 = visz + visexc;

end

function ff = FFF_Cole(ReN, Dia, eta)

RelRough = eta / Dia;

%Initial guess (multiply by 4 to give Darcy or Moody friction factor)
ff_new = 4 * F_Fact(Dia, ReN, eta, 1);
a_new = sqrt(1.0 / ff_new);

a_diff = 1E-05;  %maximum fractional difference between new  and old values for end of calculation

%Iterative calculation, counter (ic) added in case of run away calculation (breaks loop)
ic = 0;
a=0;
while ~(ic > 1000 || abs(a / a_new - 1) < a_diff)

    a = a_new;%update of calculation number

    %calculation of function
    Func_F = a + 2 * log10(RelRough / 3.7 + (2.51 / ReN) * a);

    %calculation of derivative of function
    DFunc_F = 1 + (2.18 / ReN) / (RelRough / 3.7 + (2.51 / ReN) * a);

    a_new = a - Func_F / DFunc_F; %new a
    ic = ic + 1;

end

%Note: The quantity 1/a_new^2 is the Darcy or Moody friction
%factor using the Colebrook-White equation. This
%quantity needs to be multiplied by 1/4 to give the
%Fanning friction factor, which is the desired quantity.

ff= 1.0 / (4.0 * a_new ^ 2);

end

function FFF = F_Fact(Dia, Re, eta, FF_eq)
RelRough = eta / Dia;

if FF_eq == 0

    % Equation in McCollum and Ogden (2006)
    % Note: The quantity 1/Calc_Temp^2 is the Haaland equation
    %   for the Darcy or Moody friction factor.
    %   The Fanning friction factor is 1/4
    %   the Darcy or Moody friction factor. The Fanning
    %   friction factor should be used according to
    %   McCollum and Ogden (2006).

    Calc_Temp = log10(6.91 / Re + ((eta / Dia) / 3.7) ^ 1.11);
    FFF = 1.0/ (4 * (-1.8 * Calc_Temp) ^ 2);
else
    %Equation from McCoy and Rubin (2008)
    % Note: The quantity 1/Calc_Temp^2 is the Zigrang
    %   and Sylvester equation for the Darcy or Moody friction.
    %   The Fanning friction
    %   factor is 1/4 the Darcy or Moody friction factor. The
    %   Fanning friction factor should be used according to
    %   McCoy and Rubin (2008).

    Calc_Temp = -2 * log10(RelRough / 3.7 - (5.02 / Re) * log10(RelRough / 3.7 - (5.02 / Re) * log10(RelRough / 3.7 + 13 / Re)));
    FFF = 1.0/ (4 * Calc_Temp ^ 2);
end

end


%finds maximum pipe segment length in meters
function l_max = L_seg_max(Dia,tonnes,lengthMi,capFact) 
    g = 9.80665; %acceleration due to gravity (m/s2)

    lengthM = lengthMi/ 0.621371 * 1000;
    elevation = 0; % ASSUMING 0 ELEVATION CHANGE
    etaIn = 0.0018; %Pipe inside surface roughness in inches
    etaMe = etaIn * 0.0254; %converting to meters

    p_in = 2200 * 6894.76;%converting from psig to pa
    p_out = 1200 * 6894.76;%converting from psig to pa

    q_max = (tonnes/capFact)*1000/(365*24*60*60); %maximum CO2 flow rate converted to kg/s
    temp = (53 - 32) * 5.0/9.0 + 273.15; %converting temp in to Kelvin

    p_liq_av = 0.5 * (p_out + p_in);
    Dens_CO2 = denDuanCO2(5E-05, p_liq_av * 1E-06, temp);
    Visc_CO2 = visFWVCO2(Dens_CO2, temp) * 1E-06;
    Re =4*q_max/(pi*Visc_CO2*Dia);
    %calculating fanning friction
    ff = FFF_Cole(Re,Dia, etaMe);

    a1 = 32 * ff * q_max^2/(pi^2 * Dens_CO2 * Dia^5);
    b1 = p_in - p_out;
    c1 = -g * Dens_CO2 * elevation / lengthM;

    if (a1-c1) < 0.001
        l_max = 9000000;
    else
        l_max = b1 / (a1 - c1);
    end
    if l_max > 9000000
    else
        l_max = 9000000;
    end
end


%{
function Dnl = Dia_next_low(Dia)

 if Dia <= 4 
        % There is no standard pipe with a diameter smaller than 4 inches in this model
        Dia_nom_next_low = 0;
        %Dia_in_next_low = 0;
 elseif Dia <= 6
        Dia_nom_next_low = 4;
       % Dia_in_next_low = 4;
 elseif Dia <= 8
        Dia_nom_next_low = 6;
       % Dia_in_next_low = 6;
 elseif Dia <= 10
        Dia_nom_next_low = 8;
%Dia_in_next_low = 8;
 elseif Dia <= 12
        Dia_nom_next_low = 10;
       % Dia_in_next_low = 10;
 elseif Dia <= 16
        Dia_nom_next_low = 12;
       % Dia_in_next_low = 12;
 elseif Dia <= 20
        Dia_nom_next_low = 16;
       % Dia_in_next_low = 15.162;
 elseif Dia <= 24
        Dia_nom_next_low = 20;
        %Dia_in_next_low = 18.952;
 elseif Dia <= 30
        Dia_nom_next_low = 24;
        %Dia_in_next_low = 22.742;
 elseif Dia <= 36
        Dia_nom_next_low = 30;
       % Dia_in_next_low = 28.428;
 elseif Dia <= 42
        Dia_nom_next_low = 36;
       % Dia_in_next_low = 34.114;
 elseif Dia <= 48
        Dia_nom_next_low = 42;
        %Dia_in_next_low = 39.8;
 else
        Dia_nom_next_low = 2000;
        %Dia_in_next_low = 99.9;
 end

   
 Dnl = Dia_nom_next_low;
    
    

end

%}