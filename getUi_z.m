function Ui = getUi_z(gCoeffN,N,kappa,Ntilde,Nh,R,z)

greal = gCoeffN(1:Ntilde);
gimag = gCoeffN(Ntilde+1:end);
h = 2*R/Nh;
[x,y] = meshgrid(-(Nh/2)*h:h:(Nh/2-1)*h,-(Nh/2)*h:h:(Nh/2-1)*h);
x = x-z(1)/kappa;
y = y-z(2)/kappa;
pin = 2*pi/N;
tN = (0:(N-1))*pin;
cs = cos(tN);
sn = sin(tN);
[a,b] = size(x);
currVec = zeros([size(x),Ntilde]);
repx = repmat(x,1,1,N);
repy = repmat(y,1,1,N);
repcs = zeros(1,1,N);
repsn = zeros(1,1,N);
repcs(1,1,:) = cs;
repcs = repmat(repcs,a,b,1);
repsn(1,1,:) = sn;
repsn = repmat(repsn,a,b,1);
reptN(1,1,:) = tN;
reptN = repmat(reptN,a,b,1);
int = (repx.*repcs +repy.*repsn);
for ka = -Ntilde/2 : Ntilde/2 -1
    currVec(:,:,ka+Ntilde/2+1) = sum(exp(1i * (kappa * int + ka * reptN)),3);
end
% 
gstretch = repmat(reshape(greal,1,1,length(greal)),a,b,1) + 1i*repmat(reshape(gimag,1,1,length(gimag)),a,b,1);
Ui = 1/(sqrt(2*pi))*2*pi/N * sum(currVec .* gstretch,3);
end