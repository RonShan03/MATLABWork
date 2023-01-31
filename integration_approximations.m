clear all;
format long g;
%Problem 1: 
g=@(x) exp(-x.^2);
disp("Trapezoidal Approximation with n=2 subintervals: ");
disp(trap_rule(-1, 1, 2, g));  
n_values=[2, 4, 8, 16, 32, 64, 128, 256].';
true_val = integral(g,-1,1);
disp('Table of trapezoidal approximation with all given values of n: ');
syms x
g_sym=sym(g);
g_prime_sym=diff(diff(g_sym, x));
g_prime=matlabFunction(g_prime_sym);
%Loop through and create table
disp("                         n          T(h)                      Error                    Error Ratio")
error=[];
for i=1:length(n_values)
    error(i)=abs(((2*exp(-1))*((2/n_values(i))^2)*2)/12);
if i>1
disp([n_values(i), trap_rule(-1, 1, n_values(i), g), error(i), error(i)/error(i-1)]);
else
    disp([n_values(i), trap_rule(-1, 1, n_values(i), g), error(i)]);
end
end
disp("Simpson's table--------------------------------------------------------")
%Loop through and create table
disp("                         n          S(h)                      Error                    Error Ratio")
error=[];
for i=1:length(n_values)
    error(i)=abs((12/180)*((2/n_values(i))^4)*2);
if i>1
disp([n_values(i), simp_rule(-1, 1, n_values(i), g), error(i), error(i)/error(i-1)]);
else
    disp([n_values(i), simp_rule(-1, 1, n_values(i), g), error(i)]);
end
end

disp("Trapezoid Rule with Richardson's Extrapolation table--------------------------------------------------------")
%Loop through and create table
disp("                         n          R(h)")
error=[];
for i=1:length(n_values)
    disp([n_values(i), Rich(-1, 1, n_values(i), g)]);
end

disp("More steps of Richardson's Extrapolation: ")
N_values=[2, 4, 8, 16, 32].';
disp("      N         R0()                     R1()                      R2()                      R3()  ")
for i=1:length(N_values)
     disp([N_values(i),Rich(-1,1,N_values(i),g),R1(-1,1,N_values(i),g),R2(-1,1,N_values(i),g),R3(-1,1,N_values(i),g)]);
end

%Problem 3: Approximate new Integrals: 

%Part a.)
func=@(x) sqrt(x);
func2=@(x)sin(cos(x));
error=[];
disp("Trapezoidal and Simpson's rule on sqrt(x)");
for i=1:length(n_values)
    error(i)=abs(((2*exp(-1))*((2/n_values(i))^2)*2)/12);
if i>1
disp([n_values(i), trap_rule(0, 2, n_values(i), func), error(i), error(i)/error(i-1)]);
disp([n_values(i), simp_rule(0, 2, n_values(i), func), error(i), error(i)/error(i-1)]);

else
    disp([n_values(i), trap_rule(0, 2, n_values(i), func), error(i)]);
    disp([n_values(i), simp_rule(0, 2, n_values(i), func), error(i)]);

end
end

for i=1:length(n_values)
     disp([n_values(i),Rich(-1,1,n_values(i),func),R1(-1,1,n_values(i),g),R2(-1,1,n_values(i),g),R3(-1,1,n_values(i),g)]);
end

disp("Trapezoidal and Simpson's rule on sin(cos(x))");
for i=1:length(n_values)
    error(i)=abs(((2*exp(-1))*((2/n_values(i))^2)*2)/12);
if i>1
disp([n_values(i), trap_rule(-pi/2, pi/2, n_values(i), func), error(i), error(i)/error(i-1)]);
disp([n_values(i), simp_rule(-pi/2, pi/2, n_values(i), func), error(i), error(i)/error(i-1)]);

else
    disp([n_values(i), trap_rule(-pi/2, pi/2, n_values(i), func), error(i)]);
    disp([n_values(i), simp_rule(-pi/2, pi/2, n_values(i), func), error(i)]);

end
end

%Trapezoidal Rule function: 
function [Approximation]=trap_rule(a, b, n, f)
x=linspace(a,b,n+1);
h=(b-a)/n;
Approximation=sum(h*(f(x(1:n))+f(x(2:n+1)))/2);
end

%Simpson's Rule Function: 

function[Approximation]=simp_rule(a,b,n,f)
h=b-a;
%Calculate Summation: 
xj=[];
for i=1:n;
    xj(i)=a+i*h;
end
sum=0;
for i=1:length(xj)
    if i==1 || i==length(xj)
        sum=sum+f(xj(i));
    else if i~=1 && i~=length(xj) && mod(i,2)==0
        sum=sum+4*f(xj(i));
    else
        sum=sum+2*f(xj(i));
    end
end
Approximation=(h/3)*sum;
end
end

function[Approximation]=R3(a,b,n,f)

Approximation=(1/63)*(64*R2(a/2,b/2,n,f)-R2(a,b,n,f));

end

function[Approximation]=R2(a,b,n,f)

Approximation=(1/15)*(16*R1(a/2,b/2,n,f)-R1(a,b,n,f));

end

function[Approximation]=R1(a,b,n,f)

Approximation=(1/3)*(4*Rich(a/2,b/2,n,f)-Rich(a,b,n,f));

end

%Trapezoidal Richardson Extrapolation: 

function[Approximation]=Rich(a,b,n,f)

Approximation=(1/3)*(4*trap_rule(a/2,b/2,n,f)-trap_rule(a,b,n,f));

end


