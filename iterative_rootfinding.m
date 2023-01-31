%Problem 1: 
syms x
f = @(x) atan(3*x+3);

fprintf("Problem 1: \n")
iterations(f, 0.4, 10^-10);

%Problem 2: 
f2 = @(x) sin(3*x)+1;
fprintf("Problem 2: \n")
iterations2(f2, 1, 200, 10^-10);
fprintf("This function does not iterate since the derivative at the initial" + ...
    "guess is greater than one. \n")


%Problem 3: 
f3 = @(x) sin(x)+1;
fprintf("Problem 3: \n")
iterations2(f3, 1, 200, 10^-10);
fprintf("This function actually iterates since the derivative at the initial " + ...
    "guess is less than one. This function decreases the error by a factor of " + ...
    "3 each time as opposed to a factor of 10 as the first function does.  ")



function [N, x] = iterations(g, x1, tol)
    syms x
    g_prime=diff(g(x), x);
    derivative=inline(diff(g(x),x, 'x'));
    x=[];
    if abs(derivative(x1))<1
        N=0
        while abs(x1-g(x1))>tol 
            E=abs(x1-g(x1));
            x(N+1)=x1;
            N=N+1;
            x1=g(x1);
            fprintf('Iteration %d: x Estimation=%.18f, Error=%.18f\n', N, x1, E);
        end
    else
        fprintf("Iteration method cannot be used since |g'(x1)| >= 1. Pick " + ...
            "a different initial guess. ")
    end
end


function [N, x] = iterations2(g, x1, M, tol)
    syms x
    g_prime=diff(g(x), x);
    derivative=inline(diff(g(x),x, 'x'));
    test=abs(derivative(x1));
    x=[];
    if abs(derivative(x1))<1
        N=0;
        while abs(x1-g(x1))>tol && N<=M
            E=abs(x1-g(x1));
            x(N+1)=x1;
            N=N+1;
            x1=g(x1);
            fprintf('Iteration %d: x Estimation=%.18f, Error=%.18f\n', N, x1, E);
        end
    else
        fprintf("Iteration method cannot be used since |g'(x1)| >= 1. Pick " + ...
            "a different initial guess. ")
    end
end






