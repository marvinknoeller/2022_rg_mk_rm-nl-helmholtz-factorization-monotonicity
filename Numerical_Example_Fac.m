close all
clear all
% parpool(16)

% Degree of L^2(S^1) functions
% Ntilde number of Fouriercoefficients
global Ntilde;  
Ntilde = 2^4;
% deltaStar from Banach fixpoint theorem
deltaStar = 1;
% Initial guess
global gCoeffN;
gCoeffN = zeros(2*Ntilde,1);
% number of discretization points for quadrature, also: points for FF obs.
global N;
N = 2^8;
global kappa;
kappa = 1;
tt = linspace(0,2*pi - 2*pi/Ntilde,Ntilde);
testx = [-5:.25:5];
testy = [-5:.25:5];
hold on
[XGrid,YGrid] = meshgrid(testx,testy);
surf(XGrid,YGrid,zeros(size(XGrid)),'FaceColor','none');
nt = size(XGrid,1);
XGridVec = reshape(XGrid,1,nt*nt);
YGridVec = reshape(YGrid,1,nt*nt);
zBVec = [XGridVec;YGridVec];
Nloc = N;
kappaloc = kappa;
Ntildeloc = Ntilde;
Nh = 2^7;
R = 5*kappaloc;
gCoeffNloc = gCoeffN;
%%
IniValsFin = inf(size(XGrid));
GMatFind = nan(nt,nt,2*Ntilde);
numShifts = 11;
ShiftSpace = linspace(-5,5,numShifts);
[zzVec1,zzVec2] = meshgrid(ShiftSpace,ShiftSpace);
XGridVec1 = reshape(zzVec1,1,numShifts^2);
YGridVec1 = reshape(zzVec2,1,numShifts^2);
zz = [XGridVec1;YGridVec1];
GMatCoeffFin = nan(size(XGrid,1),size(XGrid,2),2*Ntilde);
zForGMatFin = nan(size(XGrid,1),size(XGrid,2),2);
%% Scan
for zIterate = 1 : numShifts^2
    z = zz(:,zIterate);
parfor ii = 1:2*Ntilde
    gCoeffii = zeros(2*Ntilde,1);
    gCoeffii(ii) = deltaStar/sqrt(1);
    gCoeffiiMat(ii,:) = gCoeffii;
    greal = gCoeffii(1:Ntilde);
    gimag = gCoeffii(Ntilde+1:end);
    g = evaluategfun_z(greal + 1i*gimag,N,Ntilde,z); 
    UiVec(:,:,ii) = getUi_z(gCoeffii,N,kappa,Ntilde,Nh,R,z);
    gCoeffiiVec(:,ii) = gCoeffii;
    gVec(:,ii) = g;
    Fg(:,ii) =NLHH(kappa,N,@nonlinear_qh2_scaled,UiVec(:,:,ii),Nh,R,0);
    Counter(ii) = abs(2*pi/N * Fg(:,ii).'*conj(gVec(:,ii)));
end
pin = 2*pi/N;
tN = (0:(N-1))*pin;
cs = cos(tN);
sn = sin(tN);
for ii = 1:2*Ntildeloc
    parfor jj = 1:nt^2 
        zz = zBVec(:,jj);
        phiz = exp(-1i*kappaloc* (zz(1) * cs + zz(2) * sn)).';
        Denominator(jj,ii) = abs(2*pi/Nloc * gVec(:,ii).'*conj(phiz))^2;
    end
end
% Filter zeros of Denominator and set them to nan
DenoIsZero = abs(Denominator)<1e-15;
Denominator(DenoIsZero) = nan;
DenoTens = reshape(Denominator,nt,nt,[]);
for ii = 1:2*Ntilde
    ValTens(:,:,ii) = Counter(ii)./(DenoTens(:,:,ii));
end
[IniVals,minInd] = min(ValTens,[],3);
[row,col] = ndgrid(1:size(ValTens,1), 1:size(ValTens,2));
dm = ValTens(sub2ind(size(ValTens), row, col, minInd));
GMatCoeff = zeros(nt,nt,2*Ntilde);
for iG = 1:size(XGrid,1)
    for jG = 1:size(XGrid,2)
        GMatCoeff(iG,jG,:) = gCoeffiiMat(minInd(iG,jG),:);
        zForGMat(iG,jG,:) = z;
    end
end
ind = sub2ind(size(ValTens), row, col, minInd);
minIndFin = IniVals <= IniValsFin; 
minIndFinTens = repmat(minIndFin,1,1,2*Ntilde);
minIndFinTensZ = repmat(minIndFin,1,1,2);
GMatCoeffFin(minIndFinTens) = GMatCoeff(minIndFinTens);
zForGMatFin(minIndFinTensZ) = zForGMat(minIndFinTensZ);
IniValsFin = min(IniVals,IniValsFin);
end
%%
global thresh
thresh = 0.1;
options = optimoptions('fmincon','Algorithm','interior-point','MaxIterations',100,...
    'UseParallel',true,'MaxFunctionEvaluations',4000000,'ObjectiveLimit',0);
options.OutputFcn = @custom_stop_fun;
A = [];
b = [];
Aeq = [];
beq = [];
lb = -deltaStar* ones(1,2*Ntilde);
ub = deltaStar * ones(1,2*Ntilde);
Znew = IniValsFin;
percount = 1;
for iz = 1:size(XGrid,1)
    for jz = 1:size(XGrid,2)        
        if IniValsFin(iz,jz)>=thresh
            percount/size(find(IniValsFin>=thresh),1)
            percount = percount + 1;
            zz = [XGrid(iz,jz);YGrid(iz,jz)];
            zzForG = squeeze(zForGMatFin(iz,jz,:));
            gCoeffNloc = squeeze(GMatCoeffFin(iz,jz,:));       
            [gCoeffMin(iz,jz,:),val(iz,jz)] = fmincon(@(gcoeff)funhandle_zAbs(gcoeff,Nloc,kappaloc,Ntildeloc,zz,Nh,R,zzForG),gCoeffNloc,...
                A,b,Aeq,beq,lb,ub,@(gcoeff)mycon(gcoeff,deltaStar),options);
            Znew(iz,jz) = min(val(iz,jz),Znew(iz,jz));
        end              
    end
end
save("FinABS")
strforprint = "exa_FM";
Finalplots

function stop = custom_stop_fun(~, optimValues, ~)
global thresh;
if optimValues.fval <= thresh
    stop = true;
else
    stop = false;
end
end
