close all;
clear all;

% E Mo Ka = 17.44 keV

hbar  = 4.135667662e-15; % eV*s

wavelength = 0.7107488; % Mo alpha
c          = 299792458;            % m/s
nu         = c./(wavelength*10);   % 1/s
omega1     = nu*2*pi;              % radians/s
Na         = 6.02214076e23;

phi =  [1.5,4.5,5.8,6.7]*180/pi*1e-3;
tth = 16.5;
beta = tth - phi;

d   = [100 , 10, 100, 1e8]; % A
N   = {'Air','Co','Au','Si'};
Z   = [1,27,79,14];
M   = [1,58.933,196.97,28.0855];
rho = [0,8.9,19.3,2.32998];

% Roughness A
rough = [0,0,0,0];


%Total absorption cross section for MoKa radiation
%includes Rayleigh, photoelectron Compton
sigma_a     = [0 1 3.10e2 3.65e4]*1e-24; % MoKa  barn/atom 1e-24 cm^2
mu_a        = sigma_a./M*Na; % Should be for each element in layer j, so it is a matrix jxn where n is the number of different atoms in the compound
mu_photoabsorption = [0,39.53,108.7,6.068]; % [cm^2/g] Henke x-ray properties of an element

mu_lambda   = rho.*mu_photoabsorption; 
mu_lambda_a = rho.*mu_a;

mu  = ones(size(d));

% Depth resolution substrate in A
res = wavelength/10;
res_sub = wavelength/10*10*10;

[ztot,za] = generate_zmesh(res,res_sub,d);

 [R,T,A,Ep_tot,Em_tot,I,dP,Nscatt,psi,I0] = parratt_optical_xray_imd(phi,omega1,d,mu,N,'s',wavelength,Z,M,rho,mu_lambda,mu_a,true,beta,res_sub,rough,za);

markers = {'k-','k.','k--','k-.'};
markers_a = {'r-','r.','r--','r-.'};

fig = figure('Color','w');
hold on;
for i = 1:size(I,1)
  plot(ztot/10,I(i,:),markers{i},'LineWidth',1.5);
end
axis([-10 20 0 inf]);
hold on;
plot([0 0 ],[0 4],'-k','LineWidth',2);
plot([1 1 ],[0 4],'-k','LineWidth',2);
plot([11 11 ],[0 4],'-k','LineWidth',2);

xlabel('depth [nm]');
ylabel('intensity');
box off;
set(gca,'FontSize',14);


