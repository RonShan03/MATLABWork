clear all;
close all;
%Problem 1: 
% t = [0, 0.99]
% dy/dt = y^2 ; y(0) = 1;
%Approximate y(0.99) using: 

%i.) Euler's Method (n=99 steps, h=0.01), y(n+1) = y(n) + h*f(x(n), y(n))
%From IVT, 
% y = @(t) 1/(1-t);
X=linspace(0,0.99, 99);
Y=[];
yexact = 1./(1-X);    
for i=1:(length(X))
    if i==1
        Y(i)=1;
    else
        Y(i)=Y(i-1)+0.01*(Y(i-1))^2;
    end
end

figure();
hold on

plot(X,yexact);
plot(X,Y);

hold off
legend('Exact',"Euler's Approximate");

y_euler=Y(99);
ans_msg="Euler of y(0.99): %0.4f \n";
fprintf(ans_msg, y_euler);

err1=max(abs(yexact-Y));



%ii.) second order Runge Kutta Method, n=99 steps, h=0.01
y0 = 1;                  
h=0.01;                   
t = 0:h:0.99;              
yexact = 1./(1-t);     
ystar = zeros(size(t));  

ystar(1) = y0;           
for i=1:(length(t)-1)
  k1 = (ystar(i))^2;              
  y1 = ystar(i)+k1*h/2;         
  k2 = (y1)^2;
  ystar(i+1) = ystar(i) + k2*h;
end
hold on 
plot(t,yexact,t,ystar);
hold off
% legend('Exact','RK-2 Approximate');
err2=max(abs(yexact-ystar));


%iii.) fourth order Runge Kutta Method, n=99 steps, h=0.01
y0 = 1;                  
h=0.01;                   
t2 = 0:h:0.99;              
yexact = 1./(1-t2);     
ystar2 = zeros(size(t2));  

ystar2(1) = y0;           
for i=1:(length(t2)-1)
  
  k1 = (ystar2(i))^2; 
  y1 = ystar2(i)+k1*h/2;      
  
  k2 = y1^2;        
  y2 = ystar2(i)+k2*h/2;     
  
  k3 = y2^2;        
  y3 = ystar2(i)+k3*h;       
  
  k4 = y3^2;        
  
  ystar2(i+1) = ystar2(i) + (k1+2*k2+2*k3+k4)*h/6; 
end

err3=max(abs(yexact-ystar2));

hold on
plot(t2,yexact,t2,ystar2);
hold off
% legend('Exact','RK-4 Approximate');
legend('Exact', "Euler's", "RK-2", "RK-4");

ans_msg2="Max error from Euler's method: %0.4f \n";
ans_msg3="Max error from Runge Kutta 2nd Order: %0.4f \n";
ans_msg4="Max error from Runge Kutta 4th Order: %0.4f \n";

fprintf(ans_msg2, err1);
fprintf(ans_msg3, err2);
fprintf(ans_msg4, err3);
fprintf("Clearly, the Euler's method does the worst, while the 2nd order Runge Kutta does intermediately well \n and the 4th order Runge Kutta does the best by far, in terms of having the smallest error ")


%PROBLEM #2: 
n_list=[2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
ans_list_euler={};
for i=1:length(n_list)
    ans_list_euler{i}=Euler(0,(3*pi/2),1/2,i);
end
disp(ans_list_euler);

ans_list_rk4={};
for i=1:length(n_list);
    ans_list_rk4{i}=rk4(0,(3*pi/2),1/2,i);
end
disp(ans_list_rk4);


ans_list_rk2={};
for i=1:length(n_list);
    ans_list_rk2{i}=rk2(0,(3*pi/2),1/2,i);
end
disp(ans_list_rk2);

disp("n             Euler's                 RK4                  ")
 


function [Y]=Euler(t1,t2,y0,n)
h=(t2-t1)/n;
t=t1:n:t2;
Y=[];
for i=1:length(t)
    if i==1
        Y(i)=y0;
    else
        Y(i)=Y(i-1)+h*(-Y(i-1)*cos(t(i-1)));
    end
end
end



function [Y]=rk4(t1,t2,y0,n)
                 
h=(t2-t1)/n;                 
t = t1:h:t2;              
Y=zeros(size(t2));

for i=1:length(t)
  
  k1 = (Y(i))^2; 
  y1 = Y(i)+k1*h/2;      
  
  k2 = y1^2;        
  y2 = Y(i)+k2*h/2;     
  
  k3 = y2^2;        
  y3 = Y(i)+k3*h;       
  
  k4 = y3^2;        
  
  Y(i+1) = Y(i) + (k1+2*k2+2*k3+k4)*h/6; 
end


end

function [Y]=rk2(t1,t2,y0,n)
                 
h=(t2-t1)/n;                 
t = t1:h:t2;              
Y=zeros(size(t2));

for i=1:length(t)
  
  k1 = (Y(i))^2; 
  y1 = Y(i)+k1*h/2;      
  
  k2 = y1^2;        
  y2 = Y(i)+k2*h/2;     
  Y(i+1) = Y(i) + k2*h;
end


end



