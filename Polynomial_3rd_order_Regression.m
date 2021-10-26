%Polynomial 3rd order regression for temperature vs depth data
%Mohammad Rheza Zamani

clear all;
clc;
T = [21.75 ; 22.68; 25.62; 30.87; 40.5; 48.72; 63.75; 96]; %Data X 

%Kernel Matrix
z = [5^3 25 5 1; 8^3 64 8 1; 14^3 196 14 1; 21^3 441 21 1; 30^3 900 30 1; 36^3 36^2 36 1; 45^3 45^2 45 1; 60^3 60^2 60 1]; %matrix kernel

% Inversion
m=inv(z'*z)*z'*T;

%ERMS
Tcal=z*m; 
t=T-Tcal ;
misfit=sqrt(mean(t.^2));

hold on 
for i=1:60
    zreg(i)=i;
    Treg(i)=(zreg(i))^3*m(1) + (zreg(i))^2*m(2)+zreg(i)*m(3)+ m(4);
end
plot(z(:,3),T,'go') 
plot(zreg,Treg,'r'); 
plot(z(:,3),Tcal,'bo') 
xlabel('T (Celcius Degree)','Fontweight','bold')
ylabel('Z (meter)','Fontweight','bold')
legend('Temperature Data','Regression Line','Calculated Temperature','Location','Northwest')
title('Polynomial 3rd order regression for temperature vs depth Data')
grid on
hold off
