function out = A(ax,bx,sx,T)

out = exp( ( sx^2/(2*bx^2) - ax/bx )*T + (ax/bx^2 - sx^2/bx^3)*(1-exp(bx*T)) + (sx^2/(4*bx^3))*(1-exp(-2*bx*T)));

end