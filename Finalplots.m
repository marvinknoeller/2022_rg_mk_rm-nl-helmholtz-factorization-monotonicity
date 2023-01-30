newVals = IniValsFin;
surf(XGrid,YGrid,newVals,'Edgecolor','k','EdgeAlpha',.2)
f = gcf;
f.Position = [560 200 560 560];
colormap(1-gray)
view(2)
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 22;
ax.XTick = -5:5;
ax.YTick = -5:5;
ax.TickLength = [0.02, 0.04];
axis([-5,5,-5,5])
axis square
colorbar
%%
hold on
tt = linspace(0,2*pi,300);
curve_x = 1*(cos(tt) + .65*cos(2*tt)-0.65);
curve_y = 1*(1.5*sin(tt));
rotang = pi/4;
RotationMat = [cos(rotang) -sin(rotang); sin(rotang) cos(rotang)];
Res = RotationMat*[curve_x;curve_y];
curve_x = Res(1,:) + 1/2;
curve_y = Res(2,:) + 1/2;
plot3(curve_x,curve_y,3*ones(length(curve_x)),'r--','LineWidth',1.5)
grid off
% stringprint = strcat('plots/','IniVals',strforprint,"1");
% print(gcf,'-depsc',stringprint)

%%
figure
newVals = Znew;
surf(XGrid,YGrid,newVals,'Edgecolor','k','EdgeAlpha',.2)
f = gcf;
f.Position = [560 200 560 560];
colormap(1-gray)
% map = colormap;
view(2)
% colorbar
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 22;
ax.XTick = -5:5;
ax.YTick = -5:5;
ax.TickLength = [0.02, 0.04];
axis([-5,5,-5,5])
axis square
colorbar
%%
hold on
tt = linspace(0,2*pi,300);
curve_x = 1*(cos(tt) + .65*cos(2*tt)-0.65);
curve_y = 1*(1.5*sin(tt));
rotang = pi/4;
RotationMat = [cos(rotang) -sin(rotang); sin(rotang) cos(rotang)];
Res = RotationMat*[curve_x;curve_y];
curve_x = Res(1,:) + 1/2;
curve_y = Res(2,:) + 1/2;
plot3(curve_x,curve_y,3*ones(length(curve_x)),'r--','LineWidth',1.5)
grid off
% stringprint = strcat('plots/','FinalVals',strforprint,"1");
% print(gcf,'-depsc',stringprint)