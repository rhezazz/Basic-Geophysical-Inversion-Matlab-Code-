%Linear regression for temperature vs depth Data
%Mohammad Rheza Zamani
clear all;
clc;
d = importdata('data1.txt');
T = d(:,1); %Temperature Data
Z = d(:,2); % Depth Data
k = zeros(length(Z), 1);
for i = 1 : length(Z) 
    k(:, 1 ) = zeros(length(Z), 1) + 1;
end

%Kernel Matrix
z = [Z k]; 

% Inversion
m=inv(z'*z)*z'*T;

%RSME
Tcal=z*m; % Y hitung
t=T-Tcal ;
rmse=sqrt(mean(t.^2));

plot(z(:,1),T,'bo') 


hold on 
%Plotting
for i=1:80
    zreg(i)=i;
    Treg(i)=zreg(i)*m(1)+m(2);
end
plot(zreg,Treg,'g');
plot(z(:,1),Tcal,'ro') 
xlabel('T (Celcius Degree)','Fontweight','bold')
ylabel('Z (meter)','Fontweight','bold')
legend('Temperature Data','Regression Line','Calculated Temperature','Location','Northwest')
title('Linear regression for temperature vs depth Data')
grid on
hold off