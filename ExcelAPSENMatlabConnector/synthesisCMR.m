function result = synthesisCMR(ParEnt, ParReal, Qinlet, yNH3inlet, yN2inlet, yH2inlet, Tinlet, Pinlet, Qsweep, yNH3sweep, yN2sweep, yH2sweep, Tsweep, Psweep)
% define geometric 
L = ParReal(1)*100; % [m] -> [cm], length of the bed
r = ParReal(2)*100; % [m] -> [cm], radius of the bed
Nz = 1000000; %, number of discretizaiton in axial/z direction
n = ParEnt(1); % number of CMRs

% define membrane
kPGPU1 = ParEnt(2); % [GPU], ammonia permeance
seleAN = ParEnt(3); % selectivity of ammonia over nitrogen, *inf
seleAH = ParEnt(4); % selectivity of ammonia over hydrogen

% operating conditions
T = Tinlet; % [C]
pT = Pinlet; % [bar], total pressure
pPer = Psweep; % [bara], permeate pressure
% inlet
Q = Qinlet/n; % [sccm], inlet flow rate
yNH3 = yNH3inlet; % composition - NH3
yN2 = yN2inlet; % composition - N2
yH2 = yH2inlet; % composition - H2
% sweep
Qs = Qsweep/n; 
yNH3s = yNH3sweep;
yN2s = yN2sweep;
yH2s = yH2sweep;


function [res] = reactor(L, r, Nz, kPGPU1, seleAN, seleAH, T, pT, Q, yNH3, yN2, yH2, pPer, Qs, yNH3s, yN2s, yH2s)
% define util
M1 = 17.0305; % [g mol-1], MW ammmonia
M2 = 28.0134; % [g mol-1], MW nitrogen
M3 = 2.01588; % [g mol-1], MW hydrogen
R = 8.314; % [J mol-1 K-1]

% define catalyst
densitycat = 4.0; % [g cm-3] approximated using the density of alpha alumina
porosity = 0.5; % porosity

% derived
T = T + 273.15; % [C] -> [K]
pTPa = pT*1e5; %[bar] -> [Pa]
mcat = densitycat*porosity*(pi*r^2*L); % [g], amount of catalyst
dz = L/(Nz-1); % [cm]
QTP = Q*(T/273.15)*(1.01325/pT); %[sccm]->[ccm], flowrate at TP
u = QTP/( pi*r^2 )/60; %[ccm]->[cm/s], velocity at TP
c = pTPa/R/T/1e6; %[Pa]->[mol cm-3], overall molar conc., constant because isobaric
% In this case, mcatvol = density * porosity 
mcatvol = mcat/(pi*r^2*L); % [g cm-3], weight of catalyst per volume
Fmg = (2*pi*r*dz)/(pi*r^2*dz); % [cm-1]
% membrane properties
kp1 = kPGPU1*3.35e-14; % [GPU]->[mol cm-2 s-1 Pa-1]
kp2 = kp1/seleAN; % [mol cm-2 s-1 Pa-1]
kp3 = kp1/seleAH; % [mol cm-2 s-1 Pa-1]

% mass data
rhou = zeros(1, Nz); % [g cm-2 s-1], overall
rho1u = zeros(1, Nz); % [g cm-2 s-1], ammonia
rho2u = zeros(1, Nz); % [g cm-2 s-1], nitrogen
rho3u = zeros(1, Nz); % [g cm-2 s-1], hydrogen
rho1u(1) = c*yNH3*u*M1; %[mol cm-3] -> [g cm-2 s-1], mass flux of ammonia at the inlet
rho2u(1) = c*yN2*u*M2; %[mol cm-3] -> [g cm-2 s-1], mass flux of nitrogen at the inlet
rho3u(1) = c*yH2*u*M3; %[mol cm-3] -> [g cm-2 s-1], mass flux of hydrogen at the inlet

% pressure data
pNH3 = zeros(1, Nz); % [bar]
pN2 = zeros(1, Nz); % [bar]
pH2 = zeros(1, Nz);% [bar]
pNH3(1) = pT*yNH3; % [bar], partial pressure of ammonia, BC value at the inlet
pN2(1) = pT*yN2; % [bar], partial pressure of nitrogen, BC value at the inlet
pH2(1) = pT*yH2; % [bar], partial pressure of hydrogen, BC value at the inlet
pPerNH3 = pPer*yNH3s; % [bar], partial pressure of ammonia in the permeate, init
pPerN2 = pPer*yN2s; % [bar], partial pressure of nitrogen in the permeate
pPerH2 = pPer*yH2s; % [bar], partial pressure of nitrogen in the permeate

% permeation data 
J1 = zeros(1, Nz); % [mol cm-2 s-1]
J2 = zeros(1, Nz); % [mol cm-2 s-1]
J3 = zeros(1, Nz); % [mol cm-2 s-1]

% integral
for i = 1:Nz-1

    % get rate rNH3
    rate = getRateSehested(T-273.15, pNH3(i), pN2(i), pH2(i)); % [mmol gcat-1 h-1]
    rate = rate*mcatvol/1000/3600; % [mmol gcat-1 h-1] -> [mol cm-3 s-1]
    
    % get permeation flux
    J1(i) = kp1*( pNH3(i) - pPerNH3 )*1e5; % [mol cm-2 s-1]
    J2(i) = kp2*( pN2(i) - pPerN2 )*1e5; % [mol cm-2 s-1]
    J3(i) = kp3*( pH2(i) - pPerH2 )*1e5; % [mol cm-2 s-1]
    
    % update species mass flux
    rho1u(i + 1) = rho1u(i) + dz*rate*1*M1 - dz*Fmg*J1(i)*M1; % [g cm-2 s-1]
    rho2u(i + 1) = rho2u(i) + dz*rate*(-0.5)*M2 - dz*Fmg*J2(i)*M2; % [g cm-2 s-1]
    rho3u(i + 1) = rho3u(i) + dz*rate*(-1.5)*M3 - dz*Fmg*J3(i)*M3; % [g cm-2 s-1]
    
    % get mass fraction
    rhou(i) = rho1u(i) + rho2u(i) + rho3u(i); % [g cm-2 s-1]
    x1 = rho1u(i)./rhou(i);
    x2 = rho2u(i)./rhou(i);
    x3 = rho3u(i)./rhou(i);
    
    % get updated density
    rho = c/( x1/M1 + x2/M2 + x3/M3 ); % [g cm-3]
    
    % get density
    % update pressure
    rho1 = rho*x1; % [g cm-3], ammonia 
    rho2 = rho*x2; % [g cm-3], nitrogen
    rho3 = rho*x3; % [g cm-3], hydrogen
    cntpNH3 = rho1/M1*R*T*1e6/1e5;% [g cm-3] -> [bar]
    cntpN2 = rho2/M2*R*T*1e6/1e5;% [g cm-3] -> [bar]
    cntpH2 = rho3/M3*R*T*1e6/1e5;% [g cm-3] -> [bar]
    
    pNH3(i + 1) = cntpNH3;
    pN2(i + 1) = cntpN2;
    pH2(i + 1) = cntpH2;

end
   
rhou(Nz) = rho1u(Nz) + rho2u(Nz) + rho3u(Nz);

xNH3 = pNH3/pT*100; % [bar] -> mol %
xN2 = pN2/pT*100; % [bar] -> mol %
xH2 = pH2/pT*100; % [bar] -> mol %

% performance metrics
Am = (pi*2*r*L)/(Nz - 1); % [cm^2] control surface area for permeation
Ac = pi*r^2; % [cm^2], cross sectional area
% mass conservation check
conservation = ((rho1u(1) + rho2u(1) + rho3u(1))*Ac - (rho1u(Nz) + rho2u(Nz) + rho3u(Nz))*Ac - (sum(J1)*M1 + sum(J2)*M2 + sum(J3)*M3)*Am )/...
    ((rho1u(1) + rho2u(1) + rho3u(1))*Ac);
% conversion
conversion = (rho2u(1)*Ac - rho2u(Nz)*Ac - sum(J2)*Am*M2) / (rho2u(1)*Ac) * 100; % [%]
% recovery
recovery = (sum(J1)*Am*M1)/(sum(J1)*Am*M1 + rho1u(Nz)*Ac - rho1u(1)*Ac) * 100; % [%]
% Retentate composition
retentateComposition = [xNH3(Nz), xN2(Nz), xH2(Nz)];
% Permeate composition
xNH3per = sum(J1)/(sum(J1) + sum(J2) + sum(J3))*100; % [%]
xN2per = sum(J2)/(sum(J1) + sum(J2) + sum(J3))*100; % [%]
xH2per = sum(J3)/(sum(J1) + sum(J2) + sum(J3))*100; % [%]
permeateComposition = [xNH3per, xN2per, xH2per];
temp = round([conservation, conversion, recovery, permeateComposition, retentateComposition], 2);

% Permeate flow
NH3per = sum(J1)*Am/1000 + Qs*yNH3s*7.45E-7/1000; % [mol cm-2 s-1] * [cm2] -> [mol s-1] -> [kmol s-1] / [sccm] -> [kmol s-1]
N2per = sum(J2)*Am/1000 + Qs*yN2s*7.45E-7/1000; % [mol cm-2 s-1] * [cm2] -> [mol s-1] -> [kmol s-1]
H2per = sum(J3)*Am/1000 + Qs*yH2s*7.45E-7/1000; % [mol cm-2 s-1] * [cm2] -> [mol s-1] -> [kmol s-1]

% Retentate flow
NH3ret = rho1u(Nz)*Ac/M1/1000 ; % [g cm-2 s-1] -> [g s-1] -> [mol s-1] -> [kmol s-1]
N2ret = rho2u(Nz)*Ac/M2/1000 ; % [g cm-2 s-1] -> [g s-1] -> [mol s-1] -> [kmol s-1]
H2ret = rho3u(Nz)*Ac/M3/1000 ; % [g cm-2 s-1] -> [g s-1] -> [mol s-1] -> [kmol s-1]

res = [NH3per, N2per, H2per, NH3ret, N2ret, H2ret];
end

result = reactor(L, r, Nz, kPGPU1, seleAN, seleAH, T, pT, Q, yNH3, yN2, yH2, pPer, Qs, yNH3s, yN2s, yH2s) * n;


end