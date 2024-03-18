function [ztot,za] = generate_zmesh(res,res_sub,d)  

origres = res;
za{1} = linspace(d(1)-res,0,round(d(1)/res)); % this one includes the end points always.

thesum   = 0;
ztot = [];
res = origres;
% I am thinking that we cannot start in the air layer because then z
% will be weird once we get to the bottom, even if mu is zero.
% so j = 2
for j = 2:1:length(d)
    % Layer with thickness d(j) owns the coordinate at the bottom
    % surface. The upper interface belongs to the layer above.

    if j == length(d)
        res = res_sub;
    end
    %ztemp = d(j):-res:0; % this one may not include the end points.
    ztemp = linspace(d(j)-res,0,round(d(j)/res)); % this one includes the end points always.
    za{j} = ztemp;

    % We start at the surface where z=0 and work our way down

    % Once we have gone through the first layer, we take stock of the
    % position. We should have gone down exactly d(j) down but the z
    % coordinates are evenly spaced and may not be exactly d(j)
    ztot = [ ztot thesum + ztemp(end:-1:1)];
    thesum   = thesum + d(j);

end
ztot = [-za{1} ztot];
end
