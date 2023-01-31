g= @(x) 1/(3*x^2 + 5)^0.5;

N=[2, 5, 10];
step_size=@(n) (1/(2^n));
h=[];

for i=1:length(N)
    h(i)=step_size(i);
end

%PROBLEM 1: 
disp("PROBLEM 1: ");
forward_approx(1, g, h, 1);
backward_approx(1, g, h);
centered_approx(1, g, h, 1);
double_centered_approx(1, g, h);

%PROBLEM 2: 
disp("PROBLEM 2: ");
f=@(x) x*exp(x) + log(x+3);
N2=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
step_size2=@(n) (1/(2^n));
h2=[];
for i=1:length(N2)
    h2(i)=step_size2(i);
end
%Actual derivative: 
f_prime=@(x) x*exp(x) + exp(x) + 1/(x+1);
%Forwards estimated derivative: 
approx_list=forward_approx(1, f, h2, 0);

%Create list of True Errors & Error Ratio: 
true_error=[];
err_ratio=[];
for i=1:10
    true_error(i)=abs(f_prime(1)-approx_list(i));
    if i>1
        err_ratio(i)=true_error(i)/true_error(i-1);
    else
        err_ratio(i)="Not Applicable";
    end 
end

%set up and display forwards table: 
N2=N2.';
h2=h2.';
approx_list=approx_list.';
true_error=true_error.';
err_ratio=err_ratio.';

disp("Problem 2 Forward Approximation Table: ")
disp("     n         step    Df(f'(1))   Error    Error Ratio ")
disp([N2, h2, approx_list, true_error, err_ratio]);

%Forwards estimated derivative: 
approx_list2=centered_approx(1, f, h2, 0);

%Create list of True Errors & Error Ratio: 
true_error2=[];
err_ratio2=[];
for i=1:10
    true_error2(i)=abs(f_prime(1)-approx_list2(i));
    if i>1
        err_ratio2(i)=true_error(i)/true_error2(i-1);
    else
        err_ratio2(i)="Not Applicable";
    end 
end

%Set up and display centered estimated derivative: 
approx_list2=approx_list2.';
true_error2=true_error2.';
err_ratio2=err_ratio2.';

disp("Problem 2 Centered Approximation Table: ")
disp("     n         step    Df(f'(1))   Error    Error Ratio ")
disp([N2, h2, approx_list2, true_error2, err_ratio2]);

%PROBLEM 3: 
%New steps N3:
N3=N2(2:end);
%Calculate steps themselves: 
h3=h2(2:end);
%Calculate approximations
approx_list3=[];
approx_list3=problem3_approx(1, f, h3);
true_error3=[];
err_ratio3=[];
for i=2:9
    true_error3(i)=abs(f_prime(1)-approx_list3(i));
    if i>2
        err_ratio3(i)=true_error3(i)/true_error3(i-1);
    else
        err_ratio3(i)="Not Applicable";
    end 
end

%Set up and display table: 
approx_list3=approx_list3.';
true_error3=true_error3.';
err_ratio3=err_ratio3.';

disp("Problem 3 Special Approximation Table: ")
disp("     n         step    Df(f'(1))   Error    Error Ratio ")
disp([N3, h3, approx_list3, true_error3, err_ratio3]);

N3=N2(2:9);
%Calculate steps themselves: 
h3=h2(2:9);
%Calculate approximations
approx_list3=[];
approx_list3=problem3_approx(1, f, h3);
true_error3=[];
err_ratio3=[];
for i=2:8
    true_error3(i)=abs(f_prime(1)-approx_list3(i));
    if i>=3
        err_ratio3(i)=true_error3(i)/true_error3(i-1);
    else
        err_ratio3(i)="Not Applicable";
    end 
end

%Set up and display table: 
approx_list3=approx_list3.';
true_error3=true_error3.';
err_ratio3=err_ratio3.';

disp("Problem 3 Special Approximation Table: ")
disp("     n         step    Df(f'(1))   Error    Error Ratio ")
disp([N3, h3, approx_list3, true_error3, err_ratio3]);


%Part 3b.)
%There is supposed to be a large error difference on the last step for the last problem
%compared to the previous one because of what happens when step size h is
%used as a separate variable to approximate the derivative Df. 



%Forward Approximations: 

function [Approximations, Errors]=forward_approx(X, f, steps, display)
if display disp("Forward Approximations: "); end
%Approximate
    Approximations=[];
    for i=1:length(steps)
        Approximations(i)=(f(X+steps(i))-f(X))/steps(i);
    end

    if display
    disp(Approximations);
    end
%Derive theoretical derivative value
    syms x
    g_sym=sym(f);
    g_prime_sym=diff(g_sym, x);
    g_prime=matlabFunction(g_prime_sym);

    theoretical = g_prime(X);

%Calculate Error
    Errors=[];

    for i=1:length(Approximations)
         Errors(i)=abs((theoretical-Approximations(i)));
    end
    if display
    disp("Errors: ")
    disp(Errors);
    end
end


%Backward Approximations: 
function [Approximations, Errors]=backward_approx(X, f, steps)
disp("Backward Approximations: ");
%Approximate
    Approximations=[];
    for i=1:length(steps)
        Approximations(i)=(f(X)-f(X-steps(i)))/steps(i);
    end
    disp("Approximations: ");
    disp(Approximations);
%Derive theoretical derivative value
    syms x
    g_sym=sym(f);
    g_prime_sym=diff(g_sym, x);
    g_prime=matlabFunction(g_prime_sym);

    theoretical = g_prime(X);

%Calculate Error
    Errors=[];

    for i=1:length(Approximations)
         Errors(i)=abs((theoretical-Approximations(i)));
    end
    disp("Errors: ")
    disp(Errors);
end

%Centered Difference Approximations: 
function [Approximations, Errors]=centered_approx(X, f, steps, display)
if display
    disp("Centered Approximations: ");
end
%Approximate
    Approximations=[];
    for i=1:length(steps)
        Approximations(i)=(f(X+steps(i))-f(X-steps(i)))/(2*steps(i));
    end
    if display
        disp("Approximations: ");
        disp(Approximations);
    end

%Derive theoretical derivative value
    syms x
    g_sym=sym(f);
    g_prime_sym=diff(g_sym, x);
    g_prime=matlabFunction(g_prime_sym);

    theoretical = g_prime(X);

%Calculate Error
    Errors=[];

    for i=1:length(Approximations)
         Errors(i)=abs((theoretical-Approximations(i)));
    end
    
    if display
    disp("Errors: ")
    disp(Errors);
    end

end

%Second Order Centered Difference Approximations: 
function [Approximations, Errors]=double_centered_approx(X, f, steps)
disp("Second Order Centered Approximations: ");
%Approximate
    Approximations=[];
    for i=1:length(steps)
        first_approx=centered_approx(X+steps(i), f, steps, 0);
        second_approx=centered_approx(X-steps(i), f, steps, 0);
        Approximations(i)=(first_approx(i)-second_approx(i))/(2*steps(i));
    end
%Derive theoretical derivative value
    syms x
    g_sym=sym(f);
    g_prime_sym=diff(g_sym, x);
    g_prime=matlabFunction(g_prime_sym);

    theoretical = g_prime(X);

%Calculate Error
    Errors=[];

    for i=1:length(Approximations)
         Errors(i)=abs((theoretical-Approximations(i)));
    end
end


%Problem 3 function: 

function [Approximations, Errors]=problem3_approx(X, f, steps)
disp("Second Order Centered Approximations: ");
%Approximate
    Approximations=[];
    for i=1:length(steps)
        Approximations(i)=(-f(X+2*steps(i))+8*f(X+steps(i))-8*f(X-steps(i))+f(X-2*steps(i)))/(12*steps(i));
    end
    disp("Approximations: ")
    disp(Approximations);
%Derive theoretical derivative value
    syms x
    g_sym=sym(f);
    g_prime_sym=diff(g_sym, x);
    g_prime=matlabFunction(g_prime_sym);

    theoretical = g_prime(X);

%Calculate Error
    Errors=[];

    for i=1:length(Approximations)
         Errors(i)=abs((theoretical-Approximations(i)));
    end
    
    disp("Errors: ")
    disp(Errors);
end




