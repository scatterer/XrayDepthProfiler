close all;
clear all;

hs = 20e-3*1e10;
bw = 0.1e-3*1e10;
A0 = hs*bw;

% E Mo Ka = 17.44 keV

hbar  = 4.135667662e-15; % eV*s

wavelength = 0.7107488; % Mo alpha
c          = 299792458;            % m/s
nu         = c./(wavelength*10);   % 1/s
omega1     = nu*2*pi;              % radians/s
eps0 = 8.8541878128e-12;
mu0  = 4*pi*1e-7;
m    = 9.1093837015e-31;
e    = 1.602176634e-19;
h    = 6.62607015e-34;
NA   = 6.02214076e23;

r0   = e^2/(4*pi*eps0*m*c^2)*1e10;


d   = [100,30,420,500e-6*1e10]; % A 500e-6*1e10
N   = {'Air','Al2O3','VZr','Al2O3'}; 

c0  = [2/5 3/5];                   % Atomic fraction
M0  = [26.982 15.999];             % Atomic mass
Z0  = [13 8];                      % Atomic number
f2  = [0.520623E-01 0.616006E-02]; % f2
sigma0 = 2*r0*wavelength*f2*1e-16; % Cross section cm2


c1  = [0.33 0.67];                 % Atomic fraction
M1  = [50.9415 91.224];             % Atomic mass
Z1  = [23 40];                      % Atomic number
f2  = [0.547023 0.590390]; % f2
sigma1 = 2*r0*wavelength*f2*1e-16; % Cross section cm2


c2  = [2/5 3/5];                   % Atomic fraction
M2  = [26.982 15.999];             % Atomic mass
Z2  = [13 8];                      % Atomic number
f2  = [0.520623E-01 0.616006E-02]; % f2
sigma2 = 2*r0*wavelength*f2*1e-16; % Cross section cm2

Z   = [1,sum(c0.*Z0),sum(c1.*Z1),sum(c2.*Z2)];
M   = [1,sum(c0.*Z0),sum(c1.*M1),sum(c2.*M2)];

rho      = [0,2.0,5.9,3.97];       % Mass density g/cm3

rough = [3,8,6];

rough = [1,1,1];
% Gullikson x-ray booklet

mu_photoabsorption = NA./M.*[0,sum(c0.*sigma0),sum(c1.*sigma1),sum(c2.*sigma2)];

mu_lambda   = rho.*mu_photoabsorption; 
mu_lambda_a = mu_lambda;
mu_a = mu_lambda_a;


res_sub = 1000;

phi = 0:0.001:5;
mu = 1;
tth = phi*2;
beta = tth-phi;
figure(1);
set(gcf,'Color','w','Position',[1 500 900 500]);
hold off;
%plot(tth,y,'-ko');
hold on;

%plot(tth,F,'-b');
set(gca,'YScale','log');

do_efield = 0;
za = 0;
psi_d = 0;
material = 0;
                                              
[R,T,A,Ep_tot,Em_tot,I,dP,Nscatt,thepsi,I0] = parratt_optical_xray_imd(phi,omega1,d,mu,material,'s',wavelength,Z,M,rho,mu_lambda,mu_a,do_efield,psi_d,res_sub,rough,za);

R1 = GaussConv(tth,R,0.025);

 
%plot(tth,R,'-bo','LineWidth',2);
plot(tth,R,'-r','LineWidth',2);


