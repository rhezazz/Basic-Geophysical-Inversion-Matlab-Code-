%Non-linear problem : Hipocenter Determination  
%Mohammad Rheza Zamani
clc;
clear all;

%Eaethquake station position
x = [20 50 40 10];
y = [10 25 50 40];
z = [0.5 0.5 0.5 0.5];
%Synthetic Model
xs = 40;
ys = 30;
zs = -20;
%P wave velocity
Vp = 4;

% Calculating Synthetic Data 
for i = 1 : length(x)
    if(zs<0 & z(i)>0) %kondisi zs<0 & z>0
        ts(i)= sqrt((x(i)-xs)^2 +(y(i)-ys)^2 + (-z(i)-zs)^2)/Vp;
    elseif(zs>0 & z(i)<0)
        ts(i)= sqrt((x(i)-xs)^2 +(y(i)-ys)^2 + (z(i)+zs)^2)/Vp;
    else %kondisi saat 1. zs = z = 0 2. zs>0 & z>0 3. zs<0 & z<0
        ts(i)= sqrt((x(i)-xs)^2 +(y(i)-ys)^2 + (z(i)-zs)^2)/Vp;
    end
    %Add noise to data
    ts(i) = ts(i) + 0.1*(2*rand(1,1));
end

%Initial guess of the hypocenter
xa = 60;
ya = 40;
za = -30;
for iterasi = 1 : 100
    for i = 1 : length(x)
       %Jacobian matrix
        if(za>0 & z(i)<0) %za>0 & z(i)<0
            j(i,1) = (xa - x(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 + (-z(i)-za)^2)*Vp);
            j(i,2) = (ya - y(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 +(-z(i)-za)^2 )*Vp);
            j(i,3) = (-za - z(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 + (-z(i)-za)^2)*Vp);
        elseif(za<0 & z(i)>0) %za<0 & z(i)>0
            j(i,1) = (xa - x(i))/(sqrt(((xa-x(i))^2) +(ya-y(i))^2 + (za+z(i))^2)*Vp);
            j(i,2) = (ya - y(i))/(sqrt(((xa-x(i))^2) +(ya-y(i))^2 +(za+z(i))^2 )*Vp);
            j(i,3) = (za + z(i))/(sqrt(((xa-x(i))^2) +(ya-y(i))^2 + (za+z(i))^2)*Vp);
        elseif(za>0 & z(i)>0) %za>0 & z(i)>0
            j(i,1) = (xa - x(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 + (-z(i)+za)^2)*Vp);
            j(i,2) = (ya - y(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 +(-z(i)+za)^2 )*Vp);
            j(i,3) = (-za + z(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 + (-z(i)+za)^2)*Vp);
        elseif(za<=0 & z(i)<=0)%za<=0 & z(i)<=0
            j(i,1) = (xa - x(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 + (z(i)-za)^2)*Vp);
            j(i,2) = (ya - y(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 +(z(i)-za)^2 )*Vp);
            j(i,3) = (za - z(i))/(sqrt(((x(i)-xa)^2) +(y(i)-ya)^2 + (z(i)-za)^2)*Vp);
        end
        
    end
    
    %Calculating calculated data and synthetic data
    for i=1:length(x)
        if(za<0 & z(i)>0) %z<0 & z>0
            t_cal(i)= sqrt((x(i)-xa)^2 +(y(i)-ya)^2 + (-z(i)-za)^2)/Vp;
        elseif(za>0 & z(i)<0)
            t_cal(i)= sqrt((x(i)-xa)^2 +(y(i)-ya)^2 + (z(i)+za)^2)/Vp;
        else %1. zs = z = 0 2. zs>0 & z>0 3. zs<0 & z<0
            t_cal(i)= sqrt((x(i)-xa)^2 +(y(i)-ya)^2 + (z(i)-za))/Vp;
        end
        dt(i) = ts(i)-t_cal(i);
    end
    %inversion
    delta_m=(inv(j'*j))*j'*dt';
    xa=xa+delta_m(1);
    ya=ya+delta_m(2);
    za=za+delta_m(3);
    hipos(iterasi,:)=[xa(1) ya(1) za(1)];
    
    %RMSE
    for i=1:length(x)
        rms(iterasi)=(sqrt((sum((ts(i)-t_cal(i)).^2))/length(x)));
    end
%Plotting Model
figure(1);
hold on
plot3(hipos(:,1),hipos(:,2),hipos(:,3),'*b'); 
plot3(x,y,z,'ro')
plot3(xs,ys,zs,'ko') 
title(['Hypocenter Model || ',num2str(length(rms)),' Iteration || RMS : ',num2str(rms(iterasi))])
grid on;
xlabel('X (m)','Fontweight','bold'); 
ylabel('Y (m)','Fontweight','bold'); 
zlabel('Z (m)','Fontweight','bold');
legend('Hypocenter Calculated','Earthquake station','Hypocenter Synthetic')
view(3)
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