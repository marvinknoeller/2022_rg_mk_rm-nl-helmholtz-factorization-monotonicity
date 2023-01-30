function [uinfty] = NLHH(kappa,nobs,q,Ui,Nh,R,plotflag)
h = 2*R/Nh;
MaxNrOfCGIterations = 5000;
%% Grid
[Xh,Yh] = meshgrid(-(Nh/2)*h:h:(Nh/2-1)*h,-(Nh/2)*h:h:(Nh/2-1)*h);
[Q,Q0,Q2] = q(Xh/kappa,Yh/kappa);
Rp = R/kappa;
if plotflag 
    figure(20)
    imagesc([-Rp,Rp-h/kappa],[-Rp,Rp-h/kappa],real(1+Q0+Q2))
    hold on
    tp = linspace(0,2*pi);
    plot(Rp/2*cos(tp),Rp/2*sin(tp),'w','linewidth',2)
    plot(Rp*cos(tp),Rp*sin(tp),'w','linewidth',2)
    hold off
    set(gca,'YDir','normal');
    xlim([-Rp Rp-h/kappa]); ylim([-Rp Rp-h/kappa]);
    axis square
    colorbar
    title('Index of refraction n^2=1+q (real part)')
    drawnow
end
%%
pin = 2*pi/nobs;
cs = cos((0:(nobs-1))*pin);
sn = sin((0:(nobs-1))*pin);
c = getc(Nh,h);
c = h^2*c;
rhs = ToepPhi(c,Q0.*Ui);
Us_0 = CGSolve(rhs,Q0,c,Nh,MaxNrOfCGIterations);
w = FixpointIteration(c,Nh,MaxNrOfCGIterations,Q0,Q2,Ui,Us_0);
V = (Q0 + Q2.*abs(w + Us_0 + Ui).^2).*(w + Us_0 + Ui);
%% Evaluate far field
Xhcol = reshape(Xh,(Nh)^2,1);
Yhcol = reshape(Yh,(Nh)^2,1);
Vcol = reshape(V,(Nh)^2,1);
Vhelp = repmat(Vcol,1,nobs);
Vhelp = h^2*exp(-1i*(Xhcol*cs+Yhcol*sn)) .* Vhelp;
farfield = sum(Vhelp,1).';
uinfty = farfield;
end
function AX = ApplyA(X,Q,c)
    AX = X - ToepPhi(c,Q.*X);
end
function AstarX = ApplyAstar(X,Q,c)
    AstarX = X - conj(Q).*ToepPhi(conj(c),X);
end
function x = CGSolve(b,Q,c,N,NrOfCGIterations)
x = zeros(N);  % initial guess
r = b-ApplyA(x,Q,c);  % residual 
s = ApplyAstar(r,Q,c);  % resudial of the normal equations
d = s;
rho = norm(r,'fro');
iter = 1;
while rho>eps*norm(b,'fro') && iter<NrOfCGIterations
    Ad = ApplyA(d,Q,c);
    alpha = norm(s,'fro')^2/norm(Ad,'fro')^2;
    x = x+alpha*d;
    r = r-alpha*Ad;
    sneu = ApplyAstar(r,Q,c);
    beta = norm(sneu,'fro')^2/norm(s,'fro')^2;
    s = sneu;
    d = s+beta*d;   
    rho = norm(r,'fro');
    iter = iter+1;
end
if iter==NrOfCGIterations
   disp('CGNE did not converge!') 
end
end
% 

function Giter = ApplyG(wl,Q0,Q2,c,N,NrOfCGIterations,Ui,Us_0)
rhs = ToepPhi(c,Q2.*abs(wl+Us_0+Ui).^2.*(wl+Us_0+Ui));
Giter = CGSolve(rhs,Q0,c,N,NrOfCGIterations);
end

function x = FixpointIteration(c,N,NrOfCGIterations,Q0,Q2,Ui,Us_0)
w = zeros(size(Ui));
err = 1;
NonConvergence = 0;
numConvs = 1;
while err>1e-5
    wlp1 = ApplyG(w,Q0,Q2,c,N,NrOfCGIterations,Ui,Us_0);
    err = norm(wlp1-w,Inf)/norm(wlp1,Inf);
    if err>1e5 || numConvs >1000
        NonConvergence = 1;
        break
    end
    w = wlp1;
    numConvs = numConvs + 1;
end
if NonConvergence == 1
    x = NaN(size(w));
else
    x = w;
end
end