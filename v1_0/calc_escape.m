function [zlayer,dP,Abs] = calc_escape(z,I,d,mu_lambda,rho_num,beta,layer_idx,I0) % 

  % Electric field strength I is defined at every point z in the layer
  % the absorption A is defined for one or several exit angles and at
  % every point in z.

  % The intensity from a layer is given by considering each volume
  % element A0*dz/sin(incident angle) and how much attenuation has
  % happened at that depth z. Therefore we need to integrate through
  % the whole layer to find the entire intensity contributed from the layer

  % This includes also allowing the radiation to leave within the layer
  % with exp(-mu*z/sind(psid)).

  % Finally once the radiation has been counted up and propagated to the
  % surface of the layer, we propagate it through the rest of the stack
  % above. Since no NEW radiation is generated at this point in the layers
  % above we simply just attenuate the beam by sum_j exp(-muj*z/sind(psi_d))

  % Position of the considered field strength: zlayer(i)

  % Positions above the surface of our layer should be propagated through
  % with a simple exponential.
   eps0 = 8.8541878128e-12;
  mu0  = 4*pi*1e-7;
  c    = 299792458;
  h    = 6.62607015e-34;
  Z0   = 1/(eps0*c);

 
  dcum = cumsum(d(2:end));

  % assume layer_idx > 1
  if layer_idx == 2
      idx =  z < dcum(layer_idx-1) & z >= 0;
  else
      idx =  z < dcum(layer_idx-1) & z >= dcum(layer_idx-2);
  end

  
  Itop = I(:,z < dcum(layer_idx-1) & z >= 0);

  zlayer = z(idx);
  Ilayer = I(:,idx);


  [kk,ll] = size(I); % I is a function of incidence angle and z.
  
  Abs = zeros(kk,length(beta));
  dP = zeros(size(Ilayer));
  dP_top = 1;

  
  % loop over incidence angles
  for i = 1:length(beta) % loop over exit angles

      s = 0;
      % we need not transport through the present layer, layer_idx and not
      % through air either so 2:layer_idx-1
      for n = 2:layer_idx-1
          s = s + mu_lambda(n)*1e-8*d(n)/sind(beta(i));
      end
      F = exp(-s);
      Alayer   = exp(-mu_lambda(layer_idx)*1e-8.*zlayer/sind(beta(i)));

      for k = 1:kk
          %mu_lambda(layer_idx)*1e-8*1e10*
          dP(k,:) = mu_lambda(layer_idx)*1e-8*1e10*Ilayer(k,:)/(2*Z0); % Absorbed energy W /(2*Z0)? mu_lambda(layer_idx)*1e-8*1e10*
          %dP_top  = mu_lambda(2)*1e-8*1e10*Itop(k,1)/(2*Z0);

          %dP(k,:) = dP(k,:)/dP_top;

          Abs(k,i) = rho_num(layer_idx)*1e-24*1e30*F*1e-10*trapz(zlayer,dP(k,:).*Alayer); %
      end
  end
end
