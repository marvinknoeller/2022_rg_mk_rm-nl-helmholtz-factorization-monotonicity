function [q,q0,q2] = nonlinear_qh2_scaled(x,y,u)
tt = linspace(0,2*pi,1000);
curve_x = 1*(cos(tt) + .65*cos(2*tt)-0.65);
curve_y = 1*(1.5*sin(tt));
rotang = pi/4;
RotationMat = [cos(rotang) -sin(rotang); sin(rotang) cos(rotang)];
Res = RotationMat*[curve_x;curve_y];
curve_x = Res(1,:) + 1/2;
curve_y = Res(2,:) + 1/2;
[Ind1,Ind2] = inpolygon(x,y,curve_x,curve_y);
Index2 = Ind1 + Ind2;
q0supp = Index2;
q1supp = Index2;
n0 = 1.47;
q0 = n0^2 - 1;
q1 = 2.5*1e-22;
tau = 3 * 1e10;
q1 = tau^2 * q1;
if nargin == 3
    q0 = q0*q0supp;
    q2 = q1*q1supp;
    q = q0 + q2.*abs(u).^2;
else
    q0 = q0*q0supp;
    q2 = q1*q1supp;
    q = NaN;
end
end