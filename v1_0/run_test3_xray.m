close all;
clear all;

% E Mo Ka = 17.44 keV

hbar  = 4.135667662e-15; % eV*s

wavelength = 0.7107488; % Mo alpha
c          = 299792458;            % m/s
nu         = c./(wavelength*10);   % 1/s
omega1     = nu*2*pi;              % radians/s
Na         = 6.02214076e23;

% phi incidence angle in degrees
phi  = linspace(0,5,300)*180/pi*1e-3;
tth  = 2*phi;
beta = tth - phi;

d   = [200 , 700, 5000]; % A
N   = {'Air','Si','Au'};
Z   = [7,14,79];
M   = [14.007,28.0855,196.97];
rho = [0,2.32998,19.3]; %

% % Total absorption cross section for MoKa radiation 
% % includes Rayleigh, photoelectron Compton
sigma_a     = [0 3.10e2 3.65e4]*1e-24; % MoKa  barn/atom 1e-24 cm^2
mu_a        = sigma_a./M*Na; % Should be for each element in layer j, so it is a matrix jxn where n is the number of different atoms in the compound
mu_photoabsorption = [0,6.068,108.7]; % [cm^2/g] Henke x-ray properties of an element 

mu_lambda   = rho.*mu_photoabsorption; 
mu_lambda_a = rho.*mu_a;

mu  = ones(size(d));

% Roughness A
rough = [0,0,0];


% Depth resolution substrate in A
res     = 1;
res_sub = 1;

[ztot,za] = generate_zmesh(res,res_sub,d);

 [R,T,A,Ep_tot,Em_tot,I,dP,Nscatt,psi,I0] = parratt_optical_xray_imd(phi,omega1,d,mu,N,'s',wavelength,Z,M,rho,mu_lambda,mu_a,true,beta,res_sub,rough,za);

rho_a    = rho./M*Na; % g/cm3/g/mole*mole-1 = FU/cm^3


for i = 1:length(psi)
    [zlayer2,dP2,A2(i)] = calc_escape(ztot,I(i,:),d,mu_lambda,rho_a,beta(i),2,I0);
    [zlayer3,dP3,A3(i)] = calc_escape(ztot,I(i,:),d,mu_lambda,rho_a,beta(i),3,I0);
end

figure('Color','w');
plot(phi*pi/180,R*10,'k--','LineWidth',1.5);
hold on;
%plot(phi*pi/180,T*10,'k-.','LineWidth',1.5);
%plot(phi*pi/180,A*10,'-ok','LineWidth',1.5);
plot(phi*pi/180,A2/max(A2)*15,'-k','LineWidth',1.5);
plot(phi*pi/180,A3/max(A3)*7,'-.k','LineWidth',1.5);
%plot(phi*pi/180,A1+A2,'-r','LineWidth',1.5);

%axis([0 5e-3 0 15]);
set(gca,'YScale','lin');
xlabel('\psi [mrad]');
ylabel('intensity');
box off;
set(gca,'FontSize',14);

