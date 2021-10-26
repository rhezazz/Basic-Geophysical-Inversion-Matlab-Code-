%Regional and Residual gravity anomaly separation using  Trend Surface
%Analysis (TSA) 2rd order
%Mohammad Rheza Zamani
clear all;
clc;
%Function : P(x,y)=a+bx+cx^2+dy+ey^2+fxy 
x = [1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 ];
y = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 ];
cba = [42 32 20 10 0 52 40 30 20 10 62 55 50 37 20 72 69 80 52 30 82 70 60 50 40 91 82 70 60 50];
x1 = x.*x;
y1 = y.*y;
xy = x.*y;
for i = 1 : length(x) 
    k(:, 1 ) = zeros(length(x), 1) + 1;
end
%Kernel Matrix
G = [k x' x1' y' y1' xy'];
%Inversion
m=inv(G'*G)*G'*cba';

%Regional Anomalies
anomali_regional = G*m;

%Residual Anomalies
anomali_residual = cba' - anomali_regional;

% Gridding Data
xv = linspace(min(x), max(x), 20);
yv = linspace(min(y), max(y), 20); 
[Xm,Ym] = meshgrid(xv, yv);
Zm = griddata(x, y, anomali_regional, Xm, Ym); 
Z = griddata (x, y, cba, Xm, Ym); 
Zn = griddata(x,y,anomali_residual,Xm,Ym);

% Plotting data
figure(1)
subplot(3,2,1)
contourf(Xm, Ym, Zm,50,'LineStyle','none')
title('Regional Anomaly Contour')
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',8,'fontweight','bold')
subplot(3,2,2)
surf(Xm, Ym, Zm)
title('Regional Anomaly 3D View') 
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',8,'fontweight','bold')
subplot(3,2,3)
contourf(Xm, Ym, Z,50,'LineStyle','none')
title('Bouger Anomaly Contour')
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',8,'fontweight','bold')
subplot(3,2,4)
surf(Xm, Ym, Z)
title('Bouger Anomaly 3D View')
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',8,'fontweight','bold')
subplot(3,2,5)
contourf(Xm, Ym, Zn,50,'LineStyle','none')
title('Residual Anomaly Contour')
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',8,'fontweight','bold')
subplot(3,2,6)
surf(Xm, Ym, Zn)
title('Residual Anomaly 3D View')
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',8,'fontweight','bold')
