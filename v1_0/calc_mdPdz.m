function [dP,Nscatt] = calc_mdPdz(wavelength,Npz,Ndpz,Et,Er,z,mu_lambda,eps_dp,rho)
  %mulambda is in 1/A
  eps0 = 8.8541878128e-12;
  mu0  = 4*pi*1e-7;
  c    = 299792458;
  h    = 6.62607015e-34;

  
  Z0 = sqrt(mu0/eps0);
  Z0 = 1/(eps0*c);

  if mu_lambda==0 % 0/0 is 1
      A1 = 2*pi*1/(wavelength)*abs(Et)^2; %/abs(E0)^2;
      A2 = 2*pi*1/(wavelength)*abs(Er)^2; %/abs(E0)^2;
      A3 = 2*pi*1/(wavelength)*2*conj(Et)*Er; %/abs(E0)^2;
  else
      A1 = 2*pi*eps_dp/(wavelength*mu_lambda)*abs(Et)^2; %/abs(E0)^2;
      A2 = 2*pi*eps_dp/(wavelength*mu_lambda)*abs(Er)^2; %/abs(E0)^2;
      A3 = 2*pi*eps_dp/(wavelength*mu_lambda)*2*conj(Et)*Er; %/abs(E0)^2;
  end
  

  b1 = 4*pi*Ndpz/wavelength;
  b2 = -b1;
  b3 = -4*pi*1i*Npz/wavelength;

  % [E] = V/m
  % [mu] = 1/m
  % [Z0] = 1/(F/m*m/s) = s/F, [F] = C/V so [Z0] = s*V/C
  % [E^2/Z0] = V^2*C/(V*s*m^2) = V*C/[m^2s] = [J]/[m^2 * s]

  % [P] = W/m^2
  % [dP/dZ] = W/m^2/dz = W*mu = W/m

  % By dividing by the energy of a photon we get
  % photons absorbed per meter


  E_photon = h*c/(wavelength*1e-10); % [J/photon]
  if A2 == 0 % we end up with limit problems 0*inf
      dP = 1/(2*Z0)*mu_lambda*1e10*real( ... %mu_lambda*1e10 mu_lambda
          A1*exp(-b1*z)  + A3*exp(-b3*z) ); %abs(E0)^2/(2*Z0)*
  else

      dP = 1/(2*Z0)*mu_lambda*1e10*real( ... %mu_lambda*1e10 %mu_lambda
          A1*exp(-b1*z) + A2*exp(-b2*z) + A3*exp(-b3*z) ); %abs(E0)^2/(2*Z0)*
  end

 
  %dP_approx = 1/(2*Z0*h*c)*2*pi*eps_dp/(mu_lambda*1e10)*abs(Et+Er).^2;
  
  % wavelength/(h*c)*dP is equal to I0 in photons/s/m^2 * mu_lambda*real(sum(...))
  % we have I0 outside of the function so we need only do the
  % wavelength/(h*c) part.
%   Attenuator_factor = 9419;
%   Area_of_beam = 0.1e-3*15e-3;  % Could choose the slit size, since all photons should count
%   I0           = 3200*Attenuator_factor/Area_of_beam; % Flux
  
  Nscatt    = 1*dP*rho*1e24; % int(Flux/m^2*1/m * S1 ) = photons/s absorbed in the material
  % 



  % dP = mu_lambda*1e10/(2*Z0)/E_photon
end