function g = evaluategfun_z(gCoeffN,N,Ntilde,z)
pin = 2*pi/N;
tN = (0:(N-1))*pin;
g = zeros(N,1);
for ka = -Ntilde/2 : Ntilde/2 -1
    g = g + 1/sqrt(2*pi) * gCoeffN(ka+Ntilde/2+1)*exp(1i*ka*tN).';
end
g = g .* exp(-1i * (z(1) * cos(tN) + z(2) * sin(tN))).';
end