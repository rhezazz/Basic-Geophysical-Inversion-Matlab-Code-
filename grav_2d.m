%Mohammad Rheza Zamani
%2-D Gravity Linear Inversion 
clear all;
clc;
G = 6.674 * 10^-11;
conts = G*(4/3)*pi;
x = [100 300 650 950];
z = [150 200 100 200];
R = 100;
rho = [2000 9000 2000 5000];
x_titik = 0:20:1000;
%Calculating Synthetic Data
for i = 1 : length(x_titik) %
    for k = 1 : 4 
        grav(i,k) = (conts*R.^3.*z(k)/(((x_titik(i)-x(k)).^2+z(k).^2).^(3/2))).*10^5;
    end
end

%Inversion
d = grav*rho';
m = inv(grav'*grav)*grav'*d;
dcal = grav*m;

%ERMS
delta_d = d - dcal;
ERMS = sqrt(mean(delta_d.^2));


%Ploting 
subplot(2,1,1)
hold on
plot(x_titik,d,'ob')
plot(x_titik,dcal,'r')
title('2D Gravity Model || ERMS : ', ERMS)
xlabel('Distance (m)','fontweight','bold')
ylabel('Gravity Anomaly (mGal)','fontweight','bold')
legend('Forward Model','Inversion Model')
hold off
subplot(2,1,2)
plot(x,z,'ob','MarkerSize',40,'MarkerFaceColor','r')
set(gca,'Ydir','reverse')
axis([0 1000 0 250])
xlabel('X(m)','fontweight','bold')
ylabel('Z (m)','fontweight','bold')
title('Body Plotting')
