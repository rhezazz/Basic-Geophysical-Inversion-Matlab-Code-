%Polynomial 2nd order regression for temperature vs depth data
%Mohammad Rheza Zamani
clear all;
clc;
d = importdata('data2.txt');
T = d(:,1); 
Z = d(:,2); 
k = zeros(length(Z), 1);
for i = 1 : length(Z) 
    k(:, 1 ) = zeros(length(Z), 1) + 1;
end

%Kernel Matrix
z = [Z.^2 Z k];

% Inversion
m=inv(z'*z)*z'*T

%RSME
Tcal=z*m; 
t=T-Tcal ;
rmse=sqrt(mean(t.^2))

%Plotting Data
hold on 
for i=1:80
    zreg(i)=i;
    Treg(i)=(zreg(i))^2*m(1)+zreg(i)*m(2)+ m(3);
end
plot(z(:,2),T,'bo') 
plot(zreg,Treg,'g'); 
plot(z(:,2),Tcal,'ro') 
xlabel('T (Celcius Degree)','Fontweight','bold')
ylabel('Z (meter)','Fontweight','bold')
legend('Temperature Data','Regression Line','Calculated Temperature','Location','Northwest')
title('Polynomial 3rd order regression for temperature vs depth Data')
grid on
hold off