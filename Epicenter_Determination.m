%Non-linear problem : Epicenter Determination  
%Mohammad Rheza Zamani

clc;
clear all;

%Eaethquake station position
x = [20 50 40 10];
y = [10 25 50 40];

%Synthetic Model
ts = [7.1 2.8 5.0 7.9];
xs = 40;
ys = 30;

%P wave velocity
Vp = 4;

%Initial guess of the hypocenter
xa = 60;
ya = 40;
for iter = 1:100 
    for i = 1 : length(x)
         %Jacobian matrix
        j(i,1) = (xa - x(i))/(sqrt(((xa-x(i))^2) +(ya-y(i))^2)*Vp);
        j(i,2) = (ya - y(i))/(sqrt(((xa-x(i))^2) +(ya-y(i))^2)*Vp);
    end
    
    %Calculating calculated data and synthetic data
    for i = 1 : length(x)
        t_cal(i)=sqrt(((xa-x(i))^2)+((ya-y(i))^2))/Vp;
        dt(i) = ts(i)-t_cal(i);
    end
    
    %inversion
    delta_m=(inv(j'*j))*(j'*dt');
    xa=xa+delta_m(1);
    ya=ya+delta_m(2);
    epis(iter,:)=[xa(1) ya(1)];
    
    %RMSE
    for i=1:length(x)
        rms(iter)=(sqrt((sum((ts(i)-t_cal(i)).^2))/length(x)));
    end
%Plotting Hypocenter
figure(1);
hold on
plot(epis(:,1),epis(:,2),'*b'); 
plot(x,y,'ro')
plot(xs,ys,'ko') 
title(['Epicenter Model || ',num2str(length(rms)),' Iteration || RMS : ',num2str(rms(iter))])
grid on;
xlabel('X (m)','Fontweight','bold'); 
ylabel('Y (m)','Fontweight','bold'); 
legend('Epicenter Calculated','Earthquake station','Epicenter Synthetic')
xlim([0 60])
ylim([0 60])
hold off
    if rms(iter)<0.001
        break
    end
end


%Plotting RMS
figure(2);
plot(1:length(rms),rms);
xlabel('Iteration','Fontweight','bold'); 
ylabel('RMS','Fontweight','bold');
title('RMS Graphic');
grid on