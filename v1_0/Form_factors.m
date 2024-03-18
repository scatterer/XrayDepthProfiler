
function [f,f_sqrd_real, f_av_sqrd_real] = Form_factors(Q, FF, composition)

form_factor = zeros(length(composition), length(Q));
f           = zeros(length(Q), 1);
f_sqrd      = zeros(length(Q), 1);

for i = 1:length(composition)
    form_factor(i,:) = FF.c(i) + FF.f_1(i) + 1j*FF.f_2(i); 
    for j = 1:5
        gaussian_basis = FF.a(i,j).*exp( -FF.b(i,j).*(Q/(4.*pi)).^2);
        form_factor(i,:) = form_factor(i,:) + transpose(gaussian_basis);
    end
end

for i=1:length(composition)
    f      = f + transpose( composition(i).*form_factor(i,:));
    f_sqrd = f_sqrd + transpose(composition(i).*form_factor(i,:).*conj(form_factor(i,:)) );
end

f_sqrd_real    = real(f_sqrd);
f_av_sqrd_real = real(f.*conj(f));


