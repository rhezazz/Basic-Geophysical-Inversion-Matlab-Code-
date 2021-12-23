%%1-D Self Potential linear Inversion (Rumus SP forward modelling : El-Kaliboy dan Al-Gami (2009))
%Inclined sheet anomaly
%Mohammad Rheza Zamani
clear all;
clc;
%Parameter
k  = 100; %Amplitudo polarisasi
z = 15; %Kedalaman dari permukaan ke titik tengah sheet
x0 = 55;  % Jarak horizontal dari sheet
alpha = 150; %Sudut inklinasi dari sheet
a = 12; %1/2 jarak lebar dari sheet

%Jarak pengukuran
x = -100:5:100;
%Perhitungan forward modelling
for i = 1 : length(x)
    V(i) = log((((x(i)-x0)-a*cosd(alpha))^2+(z - a*sind(alpha))^2)/((((x(i)-x0)+a*cosd(alpha))^2+(z + a*sind(alpha))^2)));
end

%Inversi
d = V'*k;
m = inv(V*V')*V*d;
dcal = V'*m;

%ERMS
delta_d = d - dcal;
ERMS = sqrt(mean(delta_d.^2));

%Plot forward
figure(1)
hold on
plot(x,d,'r.','MarkerSize',20,'MarkerFaceColor','r');
plot(x,dcal,'b-')
title('1D SP Inclined Sheet Model || ERMS : ', ERMS)
legend({'Forward Model','Inversion Model'},'Location','Southeast')
xlabel('Distance (m)','FontWeight','bold');
ylabel('SP Anomaly (mV)','FontWeight','bold');
hold off
