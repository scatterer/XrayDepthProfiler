function [R,T,A,Ep_tot,Em_tot,I,dP,Nscatt,thepsi,I0] = parratt_optical_xray_imd(phi,omega,d,mu,material,sop,wavelength,Z,M,rho,mu_lambda,mu_a,do_efield,psi_d,res_sub,sigma,za)

Ep_tot = 0;
Em_tot = 0;
I = 0;
dP = 0;

Nscatt = 0;
  eps0 = 8.8541878128e-12;
  mu0  = 4*pi*1e-7;
  c    = 299792458;
  h    = 6.62607015e-34;

  Z0 = 1/(eps0*c);

thepsi = [];
% STARTING WITH ONLY TE WAVES -> S
  % kz is a vector of length 1..N 
  % d  is a vector of length 1..N thickness of each slab, d(1) is zero
  % mu is a vector of length 1..N magnetic permeability
  % phi is a vector of length 1..M is the incident grazing angle
  % omega is a vector of length 1..S is the energy in radians/sec
  
  % Henke
  % n = 1 - N*r0*lambda^2*f1/(2*pi) - N*r0*lambda^2*f2/(2*pi)
  % n = 1 - delta - 1i*beta
  % delta = N*r0*lambda^2*f1/(2*pi)
  % beta  = N*r0*lambda^2*f2/(2*pi)
  % mu    = 2*r0*lambda*f2 % cross section, not attenuation

  % de Boer 1991
  % n = epsp - 1ieps_dp
  % epsp   = 1-2*delta
  % eps_dp = 2*beta

  % epsp   = 1 - N*r0*lambda^2*f1/(pi) , where N = rho/M*Na
  % eps_dp = N*r0*lambda^2*f2/pi

  % eps_dp approx lambda*mu/(2*pi)
  % eps_dp approx n*r0*lambda^2*f2/pi= eps_dp
  % so mu_jlambda = eps_dp*2*pi/lambda
  % OK consistent.
  


  eps0 = 8.8541878128e-12;
  mu0  = 4*pi*1e-7;
  c    = 299792458;
  m    = 9.1093837015e-31;
  e    = 1.602176634e-19;
  h    = 6.62607015e-34;
  NA   = 6.02214076e23;

  r0   = e^2/(4*pi*eps0*m*c^2);

  N       = length(d) - 2; % By definition N is the number of layers,

  da = d(1);
  ds = d(end);
  d  = d(2:end-1);

  psi = zeros(size(d));
  
  
                           % not counting vacuum and substrate.
  n_phi   = length(phi);
  n_omega = length(omega);
  
 
  %
  %  Convention in the paper
  %
  %  a           vacuum       d(0)   = 0      
  %  ---------------------
  %  i = 1       layer 1      d(1)   = d1
  %  ---------------------
  %  i = i       layer i      d(i)   = di
  %  ---------------------                    z(i) = z(i-1) + d(i)
  %  j = i+1     layer i+1    d(j)   = d_(i+1)
  %  ---------------------
  %  j = N       layer N      d(N)   = dN
  %  ---------------------                   R(N)
  %  j = ns     substrate    d(ns) = 0     
  %                                          R(ns) = 0

   % n = epsp - 1ieps_dp

   na = 1;
   eps_a_p = 1;
   eps_a_dp = 0;
   mu_lambda_a = mu_lambda(1)*1e-8;
   rho_num_a = rho(1)/M(1)*NA*1e-24;

   n = zeros(N,1);
   eps_p  = zeros(N,1);
   eps_dp = zeros(N,1);

   for i = 1:N
       % real part
       eps_p(i)  = 1 - NA*r0*1e10*wavelength^2*rho(i+1)*1e-24*Z(i+1)/(pi*M(i+1));

       % imaginary part
       eps_dp(i) = wavelength*mu_lambda(i+1)*1e-8/(2*pi)*(eps_p(i) + ...
           wavelength^2*mu_lambda(i+1)^2*1e-16/(16*pi^2) )^(1/2);

       mu_lambda_indexed(i) = mu_lambda(i+1)*1e-8;
       rho_num_indexed(i)   = rho(i+1)/M(i+1)*NA*1e-24; % number of scatteres per cubic A
       n(i) = eps_p(i) - 1i*eps_dp(i);

   end

   eps_ns_p  = 1 - NA*r0*1e10*wavelength^2*rho(N+2)*1e-24*Z(N+2)/(pi*M(N+2));
   eps_ns_dp = wavelength*mu_lambda(N+2)*1e-8/(2*pi)*(eps_ns_p + ...
       wavelength^2*mu_lambda(N+2)^2*1e-16/(16*pi^2) )^(1/2);

   mu_lambda_ns = mu_lambda(end)*1e-8;
   ns = eps_ns_p - 1i*eps_ns_dp;
   rho_num_ns = rho(end)/M(end)*NA*1e-24;

   n = [n ;ns];

   eps_p = [eps_p;eps_ns_p];
   eps_dp = [eps_dp; eps_ns_dp];

   R    = zeros(n_phi,1);
   T    = zeros(n_phi,1);
   A    = zeros(n_phi,1);
  % Ang_part = [];
   I = [];
   dP = [];
   Nscatt = [];

  for k = 1:n_omega
      for l = 1:n_phi


          % %
          % r_ns_ns+1+rns+1*exp(2*1i*2*pi/wavelength*d_ns*n_ns*cosd(theta_ns))/
          % 1 + r_ns_ns+1*r_ns*phase but r_ns_ns+1  is reflection
          % coefficient at the back of the substrate

          % Let us start with an infinite substrate. Then 
          % r_N+1 = r_N+1,N+2 because r_N+2 = 0 nothing reflects back.
          % t_N+1 = 0, because t_N+2 = 0 because nothing transmits through


          r = zeros(N,1);
          t = zeros(N,1);

          r(N+1) = 0; % r(ns) is the net reflection from the bottom of the substrate
          t(N+1) = 1; % t(ns) is the net transmission from the bottom of the substrate

          thetaa   = (90-phi(l))*pi/180;
          theta    = zeros(N,1);
          phase    = zeros(N,1);
          beta     = zeros(N,1);
          Nz       = zeros(N,1);
          Nzp      = zeros(N,1);
          Nzdp     = zeros(N,1);
          psi      = zeros(N,1);

          theta(1) = asin(na/n(1)*sin(thetaa));

          Nzpa   = sqrt( 0.5*(1 - cosd(phi(l))^2 + sqrt( (1 - cosd(phi(l))^2)^2 + 0^2)) );
          Nzdpa  = 0;


          Nz(1)    = sqrt( n(1) - cosd(phi(l))^2 );
          Nzp(1)   = sqrt( 0.5*(eps_p(1) - cosd(phi(l))^2 + sqrt( (eps_p(1) - cosd(phi(l))^2)^2 + eps_dp(1)^2)) );
          Nzdp(1)  = eps_dp(1)/(2*Nzp(1));
          Nxp      = cosd(phi);
          psi(1)   = asind(Nzp(1)); % atand(Nzp(1)/Nxp);



          beta(1)  = pi*d(1)*Nz(1)/wavelength;
          phase(1) = exp(-2*1i*beta(1));
          for j = 2:N
            theta(j) = asin(n(j-1)/n(j)*sin(theta(j-1)));
            Nz(j)    = sqrt( n(j) - cosd(phi(l))^2 );
            Nzp(j)   = sqrt( 0.5*(eps_p(j) - cosd(phi(l))^2 + sqrt( (eps_p(j) - cosd(phi(l))^2)^2 + eps_dp(j)^2)) );
            Nzdp(j)  = eps_dp(j)/(2*Nzp(j));

            beta(j)  = pi*d(j)*Nz(j)/wavelength;
            phase(j) = exp(-2*1i*beta(j));
            psi(j)   = asind(Nzp(j));
          end

          thetas = asin(n(N)/ns*sin(theta(N)));
          Ns = sqrt( ns - cosd(phi(l))^2 );
          Nszp   = sqrt( 0.5*(eps_ns_p(1) - cosd(phi(l))^2 + sqrt( (eps_ns_p(1) - cosd(phi(l))^2)^2 + eps_ns_dp(1)^2)) );
          Nszdp  = eps_ns_dp(1)/(2*Nszp(1));
          psi(N+1) = asind(Nszp);

          betas  = 2*pi*ds*Ns/wavelength;
          phases = exp(-1i*betas);

          phase  = [phase; phases];
          theta  = [theta; thetas];
          Nz     = [Nz; Ns];
          Nzp    = [Nzp; Nszp];
          Nzdp   = [Nzdp; Nszdp];

          Ep_tot = [];
          Em_tot = [];
          dPz_tot_alt = [];
          scatterers = [];
          %tot_ang_part = [];
          if do_efield

              Ep_tot(1) = 1;
              Em_tot(1) = 0;
              Ep_i = 0;
              Em_i = 0;
              

              z = za{end};
              betaz = 2*pi*Ns/wavelength*z(end:-1:1);
              Ep_tot  = 1*exp(-1i*betaz);
              Em_tot  = zeros(1,length(z));

              %P1z = 1/(4*Z0)*(Ep_tot+Em_tot).*(Ep_tot-Em_tot)*Ns;
              %P2z = conj(P1z);
              %Pz_tot = P1z + P2z;
              %curdP  = Pz_tot;
              curdP    = mu_lambda_ns*1e10/(2*Z0)*abs(Ep_tot + Em_tot).^2;
              curscatt = 1*curdP*rho_num_ns*1e24; % int(Flux/m^2*1/m * S1 ) = photons/s absorbed in the material


              Et = 1;
              Er = 0;

             %[temp,curscatt] = calc_mdPdz(wavelength,Nszp,Nszdp,Et,Er,z(end:-1:1),mu_lambda_ns,eps_ns_dp,rho_num_ns);
             dPz_tot_alt = [curdP dPz_tot_alt];
             scatterers  = [curscatt scatterers];

          
             
          end

          for i = N:-1:1

              j = i + 1;
              
              %r_ij = (n(i)*cos(theta(i))-n(j)*cos(theta(j)))/(n(i)*cos(theta(i))+n(j)*cos(theta(j)))
              %t_ij = 2*n(i)*cos(theta(i))/(n(i)*cos(theta(i))+n(j)*cos(theta(j)));

              S(i)  = exp(-2*(2*pi*sigma(j)/wavelength)^2*Nz(i)*Nz(j));
              Sp(i) = exp( (2*pi*sigma(j)/wavelength)^2*(Nz(i) - Nz(j))^2/2); 
              r_ij = (Nz(i)-Nz(j))/(Nz(i)+Nz(j))*S(i);
              t_ij = 2*Nz(i)/(Nz(i)+Nz(j))*Sp(i);
     
              r(i) = (r_ij + r(j)*phase(j)^2)/(1 + r_ij*r(j)*phase(j)^2);
              t(i) =     t_ij*t(j)*phase(j)^2/(1 + r_ij*r(j)*phase(j)^2);

              if do_efield
                  z = za{i+1};
                  betaz = -2*pi*Nz(i)/wavelength*z;
                  Ep_i  =    1/t_ij*exp(-1i*betaz)*Ep_tot(1) + r_ij/t_ij*exp(-1i*betaz)*Em_tot(1);
                  Em_i  = r_ij/t_ij*exp( 1i*betaz)*Ep_tot(1) + 1/t_ij   *exp( 1i*betaz     )*Em_tot(1);

                  %P1z = 1/(4*Z0)*(Ep_i+Em_i).*(Em_i-Ep_i).*Nz(j);
                  %P2z = conj(P1z);
                  %Pz_tot = P1z + P2z;

                  %curdP  = Pz_tot; %diff(Pz_tot)./diff(z);
                  curdP = mu_lambda_indexed(i)*1e10/(2*Z0)*abs(Ep_i + Em_i).^2;
                  curscatt = 1*curdP*rho_num_indexed(i)*1e24; % int(Flux/m^2*1/m * S1 ) = photons/s absorbed in the material

                  % Same as eq 4 in PRB 38, 8579 1988. Transfer matrix
                  Er = 1/t_ij*Ep_tot(1) + r_ij/t_ij*Em_tot(1);
                  Et = r_ij/t_ij*Ep_tot(1) + 1/t_ij*Em_tot(1);

                  %Er  = phase(i)^2*r(i)*Et ?
                  %Et  = t(j)/t(i)*Et_j+1 ?
                
                 % tot_ang_part = [ang_part tot_ang_part];
                  Ep_tot = [Ep_i Ep_tot];
                  Em_tot = [Em_i Em_tot];


                 % [temp,curscatt] = calc_mdPdz(wavelength,Nzp(i),Nzdp(i),Et,Er,z,mu_lambda_indexed(i),eps_dp(i),rho_num_indexed(i));

             %        figure(1)
             % plot(temp*1e-10,'-bo')
             % hold on
             % plot(curdP,'-xr')
             % set(gca,'YScale','log')
                  
                  dPz_tot_alt = [curdP dPz_tot_alt];  
                  scatterers  = [curscatt scatterers];

              end
          end


          Na = sqrt( na - cosd(phi(l))^2 );

          %r_a1 = (na*cos(thetaa)-n(1)*cos(theta(1)))/(na*cos(thetaa)+n(1)*cos(theta(1)))
          %t_a1 = 2*na*cos(thetaa)/(na*cos(thetaa)+n(1)*cos(theta(1)));
          Sa  = exp(-2*(2*pi*sigma(1)/wavelength)^2*Na*Nz(1));
          Sap = exp( (2*pi*sigma(1)/wavelength)^2*(Na - Nz(1))^2/2); 
          
          r_a1 = (Na-Nz(1))/(Na+Nz(1))*Sa;
          t_a1 = 2*Na/(Na+Nz(1))*Sap;

          
          r_a = (r_a1 + r(1)*phase(1)^2)/(1 + r_a1*r(1)*phase(1)^2);
          t_a = t_a1*t(1)*phase(1)^2/(1 + r_a1*r(1)*phase(1)^2);
          
          if do_efield

              
              betaz = -2*pi*Na/wavelength*za{1};

              % Ep_tot(1) is the top of the first real layer, i.e. surface
              Ep_a  =    1/t_a1*exp(-1i*betaz)*Ep_tot(1) + r_a1/t_a1*exp(-1i*betaz)*Em_tot(1);
              Em_a  = r_a1/t_a1*exp( 1i*betaz)*Ep_tot(1) + 1/t_a1   *exp( 1i*betaz)*Em_tot(1);
              
              %P1z = 1/(4*Z0)*(Ep_a+Em_a).*(Ep_a-Em_a)*Na;
              %P2z = conj(P1z);
              %Pz_tot = P1z + P2z;
              %curdP  = Pz_tot; %diff(Pz_tot)./diff(za{1});
              curdP = mu_lambda_a*1e10/(2*Z0)*abs(Ep_a + Em_a).^2;
              curscatt = 1*curdP*rho_num_a*1e24; % int(Flux/m^2*1/m * S1 ) = photons/s absorbed in the material

              Ep_aa  =    1/t_a1*Ep_tot(1) + r_a1/t_a1*Em_tot(1); 
              Em_aa  = r_a1/t_a1*Ep_tot(1) + 1/t_a1   *Em_tot(1);

              Er = (1/t_a1*Ep_tot(1) + r_a1/t_a1*Em_tot(1));
              Et = (r_a1/t_a1*Ep_tot(1) + 1/t_a1*Em_tot(1));

              Ep_tot = [Ep_a Ep_tot];
              Em_tot = [Em_a Em_tot];

              % We need the Electric field at the top of the stack, not of the standing
              % wave somewhere in vacuum. So last of the air layer.
              r_a_alt = Em_aa(end)/Ep_aa(end);

              t_a_alt = 1/Ep_aa(end);

             

                
              Ep_tot = t_a_alt*Ep_tot;
              Em_tot = t_a_alt*Em_tot;

             
              I(l,:) = abs(Ep_tot+Em_tot).^2;
              I0     = t_a_alt;

               
             % [temp,curscatt] = calc_mdPdz(wavelength,Nzpa,Nzdpa,Et,Er,za{1},mu_lambda_a,eps_a_dp,rho_num_a);

              dPz_tot_alt = [curdP dPz_tot_alt]; % remember mu_lambda is not the same indexing
              scatterers  = [curscatt scatterers];

              
              dP(l,:) = abs(t_a_alt)^2*dPz_tot_alt; % Here we multiply by one divided by the electric field strength at z=0
              
              Nscatt(l,:) = abs(t_a_alt)^2*scatterers;

              thepsi(l,:) = psi;

              % The emitted flourescence is
              
              % I = lambda/(h*c) * int(dPdz)
              % But dPdZ is E0^2/(2*Z0)
              % so I = I0*sum(A*int()) instead. and the A coefficients are
              % divided by E0^2.

              % The absorption factor is the number of absorbed photons
              % divided by the number of incoming photons, per second, per
              % surface area, per solid angle
              %I0 = wavelength/(h*c); % /divide by 4pi to get per solid angle and multiply by abs(Ep_aa(end))^2 to get the electric field. 
              % 1/(2Z0) is already in P0
              % Because we are dividing by I0 and we have already accounted
              % for abs(Ep_aa(end))^2 in dP we shouldnt do it again.

              % C^a_nj - atomic fraction of n in j
              % C_aj - mass fraction of element a in j
              %
              % rho_j - mass density in g/A3

              %absorbed_photons(l,:) = dP(l,:)/I0;

              %dPz_tot_alt = [sqrt(mu_lambda_a*ones(size(za{1}))) dPz_tot_alt];
              %dP(l,:) = I(l,:).*dPz_tot_alt.^2;
              % it seems we need to build the ang_part from the surface and
              % down and not bottom up like the electric field.

             % ang_part = exp(-mu_lambda(i)*1e-8/sind(psi_d)*za{1}(end:-1:1)); % should be mu_a
                 
             % Ang_part(l,:) = [ang_part tot_ang_part];
          else
              I0 = 1;
          end


          t_s = prod([t_a; t(1:end-1)]);

          R(l,k) = abs(r_a)^2;
          T(l,k) = real(Ns/Na)*abs(t_s)^2; 

         
          T(isnan(T(l,k)),k) =0;
          A(l,k) = 1 - R(l,k) - T(l,k);
          %r  = Em(1,1)/Ep(1,1);
          %t  = 1/Ep(1,1);

          %Ep = Ep*t;
          %Em = Em*t;
      end
  end
end