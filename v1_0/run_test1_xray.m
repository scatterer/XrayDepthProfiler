close all;
clear all;

% E Mo Ka = 17.44 keV
hbar  = 4.135667662e-15; % eV*s
   eps0 = 8.8541878128e-12;
  mu0  = 4*pi*1e-7;
  c    = 299792458;
  h    = 6.62607015e-34;
  Z0   = 1/(eps0*c);
hbar  = 4.135667662e-15; % eV*s

wavelength = 0.7107488; % Mo alpha
c          = 299792458;            % m/s
nu         = c./(wavelength*10);   % 1/s
omega1     = nu*2*pi;              % radians/s
Na         = 6.02214076e23;

% phi incidence angle in degrees
phi =  [1.7,1.85,1.92,2.3]*180/pi*1e-3;

% tth is twotheta
tth = 16.5;

% beta is detector angle
beta = tth - phi;

% Thickness in A
d   = [200 , 700, 1e7]; % A
N   = {'Air','Si','Au'};
% Atomic number
Z   = [7,14,79];
% Molecular weight [g/mole]
M   = [14.007,28.0855,196.97];
% Mass density in [g/cm3]
rho = [0,2.32998,19.3]; %2.32998

% Roughness A
rough = [0,0,0];


% Total absorption cross section for MoKa radiation 
% includes Rayleigh, photoelectron Compton for the element/compound
sigma_a     = [0 3.10e2 3.65e4]*1e-24; % MoKa  barn/atom 1e-24 cm^2
mu_a        = sigma_a./M*Na; % Should be for each element in layer j, so it is a matrix jxn where n is the number of different atoms in the compound

mu_photoabsorption = [0,6.068,108.7]; % [cm^2/g] Henke x-ray properties of the compound

mu_lambda   = rho.*mu_photoabsorption; 

mu  = ones(size(d));

res = 1;
% Depth resolution substrate in A
res_sub = 10;

% only s is implemented but p is very similar for x-rays  
[ztot,za] = generate_zmesh(res,res_sub,d);
%[R,T,A,Ep_tot,Em_tot,I,Ang_part,dP,Nscatt] = parratt_optical_xray_imd(phi,omega1,d,mu,N,'s',wavelength,Z,M,rho,mu_lambda,mu_a,true,beta,res_sub,rough,za);
[R,T,A,Ep_tot,Em_tot,I,dP,Nscatt,psi,I0] = parratt_optical_xray_imd(phi,omega1,d,mu,N,'s',wavelength,Z,M,rho,mu_lambda,mu_a,true,beta,res_sub,rough,za);

%z = 1:1:size(I,2); 
%z = z -d(1);

markers = {'k-','k.','k--','k-.'};

fig = figure('Color','w');
hold on;
for i = 1:size(I,1)
  plot(ztot/10,I(i,:),markers{i},'LineWidth',1.5);

 % plot(z/10,dP(i,:),markers{i},'LineWidth',1.5,'Color','r');
end
plot([0 0 ],[0 30],'-k','LineWidth',2);
plot([70 70 ],[0 30],'-k','LineWidth',2);

axis([-20 80 0 inf]);
xlabel('depth [nm]');
ylabel('intensity');
box off;
set(gca,'FontSize',14);
