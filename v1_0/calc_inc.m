function I = calc_inc(lambda,c,Z,Q)
h      = 6.62607004e-34;
m      = 9.10938356e-31;
c0     = 299792458;

theta = asind(Q/4/pi*lambda);
lambdap = lambda + 2*h/(m*c0)*sind(theta).^2*1e10;

n = length(c);
i = 0;

a = 2.6917*Z.^(-1) + 1.2450;
b = 1.1870*Z.^(-1) + 0.1075 + 0.00436*Z - (0.01543*Z).^2 + (0.01422*Z).^3;

for j = 1:n
    i = i + c(j)*Z(j)*(b(j)*Q).^a(j)./(1+(b(j).*Q).^a(j));
end

I = (lambda./lambdap).^2.*i; %

end