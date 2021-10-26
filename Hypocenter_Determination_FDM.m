%Non-linear problem : Hipocenter Determination using Finite Difference
%Method
%Mohammad Rheza Zamani
clear all 
clc;

%Earthquake station position
xs = [10 50 10 50];
ys = [20 20 60 60];
zs = [0.5 0.5 0.5 0.5];

%Synthetic Model
xg = 40;
yg = 30;
zg = 20;
%P wave velocity
v=5;
% Calculating Synthetic Data 
for i=1: length(xs)
    t(i) = ((((xs(i)-xg)^2+(ys(i)-yg)^2+(zs(i)-zg)^2))^0.5)/v;
    t(i) = t(i)+0.1*((2*rand(1,1))-1);
end

%Initial guess of the hypocenter
x0= 30;
y0= 30;
z0= 30;
m0= [x0 y0 z0];
for iterasi = 1:100
%Jacobian matrix
%Pertubation = 10%
    for i=1:length(xs)
        t0(i) = (((xs(i)-x0)^2+(ys(i)-y0)^2+(zs(i)-z0)^2)^0.5)/v; 
        %Element x for jacobian matrix parameter model
        tx(i) = ((xs(i)-1.1*x0)^2+(ys(i)-y0)^2+(zs(i)-z0)^2)^0.5/v;       
        j(i,1)= (tx(i)-t0(i))/(0.1*x0);
        %Element y for jacobian matrix parameter model
        ty(i) = ((xs(i)-x0)^2+(ys(i)-1.1*y0)^2+(zs(i)-z0)^2)^0.5/v;  
        j(i,2)= (ty(i)-t0(i))/(0.1*y0);
        %Element z for jacobian matrix parameter model
        tz(i) = ((xs(i)-x0)^2+(ys(i)-y0)^2+(zs(i)-1.1*z0)^2)^0.5/v;          
        j(i,3)= (tz(i)-t0(i))/(0.1*z0);
    end
%Calculating calculated data and synthetic data difference
    for i=1:length(xs)
        t_cal(i)=(sqrt(((x0-xs(i))^2)+((y0-ys(i))^2)+((z0-zs(i))^2)))/v;
        dt(i) = t(i)-t_cal(i);
    end
%Inversion
    deltam=(inv(j'*j))*(j'*dt');
    x0=x0+deltam(1);
    y0=y0+deltam(2);
    z0=z0+deltam(3);
    hipos(iterasi,:)=[x0(1) y0(1) z0(1)];
    
%RMSE
    rms(iterasi)=0;
    for i=1:1:length(xs)
    rms(iterasi)=rms(iterasi)+(sqrt((dt(i)-mean(dt))^2))/(length(xs)-1);
    end
%plot hiposenter
figure(1)
hold on
plot3(hipos(:,1),hipos(:,2),hipos(:,3),'*b');
plot3(xs,ys,zs,'ro')
plot3(xg,yg,zg,'ko')
title(['Hypocenter Model || ',num2str(length(rms)),' Iteration || RMS : ',num2str(rms(iterasi))])
xlabel('X (m)','Fontweight','bold'); 
ylabel('Y (m)','Fontweight','bold'); 
zlabel('Z (m)','Fontweight','bold');
legend('Hypocenter Calculated','Earthquake station','Hypocenter Synthetic')
set(gca,'Zdir','reverse')
grid on
view(3)
hold off
    if rms(iterasi)<0.001
        break
    end
end

%Plotting rms
figure(2);
plot(1:length(rms),rms);
xlabel('Iteration','Fontweight','bold'); 
ylabel('RMS','Fontweight','bold');
title('RMS Graphic');
grid on