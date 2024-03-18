function Fpad = padme(F,thesize)

  half=thesize/2;
  Fpad                       = zeros(length(F)+half*2,1);
  Fpad(1:half)               = ones(half,1).*F(1);
  Fpad(half+1:end-half)      = F;
  Fpad(end-half+1:end)       = ones(half,1).*F(end);
end
