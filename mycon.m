function [c, ceq] = mycon(g,deltaStar)
N = length(g)/2;
ceq = sqrt(sum(abs(g(1:N)+1i*g(N+1:end)).^2))-deltaStar;

c = [];

end