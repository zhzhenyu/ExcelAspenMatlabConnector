function r = getRateSehested(T, pNH3, pN2, pH2)
% derived
T = T + 273.15; % C -> K

% kinetics parameters
R = 8.314; % [J mol-1 K-1]
Kb = 2.16E3*exp(-48*1000/R/T);
Ka = 2.73E-2*exp(-27.1*1000/R/T);
% 2NsK1k2 renamed k
k = 52*150*exp(-6.6*1000/R/T);
Keq = 2.03E-12*exp(101.6*1000/R/T);

% get rate
r = k*(pN2 - pNH3.^2./(pH2.^3*Keq)).*(1 + pNH3./(pH2.^1.5*Ka) + pH2.^0.5/Kb).^-2; %[umol g-1 s-1]
r = r*3600/1000; % [umol g-1 s-1] -> [mmol g-1 h-1]
end