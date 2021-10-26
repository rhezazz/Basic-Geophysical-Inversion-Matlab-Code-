%Mohammad Rheza Zamani
%3D Gravity Linear Inversion 
clear all;
clc;
G = 6.674 * 10^-11;
conts = G*(4/3)*pi;
x = [100 300 650 950];
y = [200 600 200 800];
z = [150 200 100 200];
R = 100;
rho = [2000 9000 2000 5000];
x_titik = 0:20:1000;
y_titik = 0:20:1000;
[x_grid y_grid] = meshgrid(x_titik,y_titik);
xi = reshape(x_grid,[51*51,1]); 
yi = reshape(y_grid,[51*51,1]);
grav = zeros(51, 4);
%Calculating Synthetic Data
for i = 1 : length(xi) 
    for k = 1 : 4 
        grav(i,k) = (conts*R.^3.*z(k)/(((xi(i)-x(k)).^2+(yi(i)-y(k)).^2+z(k).^2).^(3/2))).*10^5;
    end
end

%Inversion
d = grav*rho';
m = inv(grav'*grav)*grav'*d;
dcal = grav*m;

%ERMS
delta_d = d - dcal;
ERMS = sqrt(mean(delta_d.^2));

%Plotting 
d_shape = reshape(d,[51,51]);
d_inv_shape = reshape(dcal,[51,51]);
subplot(1,2,1)
contourf(x_grid,y_grid,d_shape,50,'LineStyle','none')
title('Synthetic Data')
xlabel('X (m)','fontweight','bold')
ylabel('Y (m)','fontweight','bold')
colormap(jet)
colorbar
title(colorbar,'mGal','FontSize',8,'fontweight','bold')
subplot(1,2,2)
contourf(x_grid,y_grid,d_inv_shape,50,'LineStyle','none')
title('Inversion Model || ERMS : ', ERMS)
xlabel('X (m)','fontweight','bold')
ylabel('Y (m)','fontweight','bold')
colorbar
title(colorbar,'mGal','FontSize',8,'fontweight','bold')
colormap(jet)


