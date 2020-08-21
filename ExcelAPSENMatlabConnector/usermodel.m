
function [ParEnt,ParReal,CorSal]=usermodel(ParEnt,ParReal,CorEnt)

matrixSize = size(CorEnt);
n = matrixSize(1) - 9; %Number of components
CorSal = CorEnt; %It makes the Output stream matrix equal to the Inlet stream matrix

% inlet info
Qinlet = CorEnt(n+1,1)*1000*1/7.45E-7; %[kmol s-1] -> [mol s-1] -> [sccm]
yNH3inlet = CorEnt(n,1)/CorEnt(n+1,1); % inlet NH3 composition
yN2inlet = CorEnt(n-2,1)/CorEnt(n+1,1); % inlet N2 composition
yH2inlet = CorEnt(n-1,1)/CorEnt(n+1,1); % inlet H2 composition
Tinlet = CorEnt(n+2,1) - 273.15; %[K] -> [C]
Pinlet = CorEnt(n+3,1)/1E5; %[N m-2] -> [bar]
% sweep info
Qsweep = CorEnt(n+1,2)*1000*1/7.45E-7; %[kmol s-1] -> [mol s-1] -> [sccm]
yNH3sweep = CorEnt(n,2)/CorEnt(n+1,2); % sweep NH3 composition
yN2sweep = CorEnt(n-2,2)/CorEnt(n+1,2); % inlet N2 composition
yH2sweep = CorEnt(n-1,2)/CorEnt(n+1,2); % inlet H2 composition
Tsweep = CorEnt(n+2,2) - 273.15; %[K] -> [C]
Psweep = CorEnt(n+3,2)/1E5; %[N m-2] -> [bar]
% call reactor model
res = synthesisCMR(ParEnt, ParReal, Qinlet, yNH3inlet, yN2inlet, yH2inlet, Tinlet, Pinlet, Qsweep, yNH3sweep, yN2sweep, yH2sweep, Tsweep, Psweep);

% write results
% ret
CorSal(n,1) = res(4); % NH3 flow in the ret
CorSal(n-2,1) = res(5); % N2 flow in the ret
CorSal(n-1,1) = res(6); % H2 flow in the ret
CorSal(n+1,1) = CorSal(n,1) + CorSal(n-1,1) + CorSal(n-2,1); % Ret flow

% per
CorSal(n,2) = res(1); % NH3 flow in the perm
CorSal(n-2,2) = res(2); % N2 flow in the perm
CorSal(n-1,2) = res(3); % H2 flow in the perm
CorSal(n+1,2) = CorSal(n,2) + CorSal(n-1,2) + CorSal(n-2,2); % Perm flow

end