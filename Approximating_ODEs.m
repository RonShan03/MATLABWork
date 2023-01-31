%Approximate y if: 
% dy/dt + 2y = 2 - e^(-4t) ==> dy/dt = -2y + 2 - e^(-4t)
% y(0) = 1
% h=0.05; t=[0,5]

f = @(t) 1 + 0.5*exp(-4*t) - 0.5*exp((-2)*t);

%h=0.05:

x=linspace(0,5,100);
y=[];
f_array=[];
error=0;
for i=1:length(x)
    if x(i)==0
        y(i)=1;
        f_array(i)=1;
    else
        y(i)=y(i-1) + 0.05*((-2)*(y(i-1)) + 2 - exp((-4)*x(i)));
        f_array(i)=f(x(i));
    end
    if abs(f_array(i)-y(i))>error

        error=abs(f_array(i)-y(i));

    end
end
error_msg='Max error: %0.4f';
fprintf(error_msg, error);

figure();
p1=plot(x,y);
title("h=0.05: ");

hold on
p2=plot(x,f_array);
hold off

legend([p1 p2],{'Approximation','Function'});

%Re-do with h=0.01: 
x=linspace(0,5,500);
figure();
y2=[];

error=0;
for i=1:length(x)
    if x(i)==0
        y2(i)=1;
        f_array(i)=1;
    else
        y2(i)=y2(i-1) + 0.01*((-2)*(y2(i-1)) + 2 - exp((-4)*x(i)));
        f_array(i)=f(x(i));
    end
    if abs(f_array(i)-y2(i))>error

        error=abs(f_array(i)-y(i));

    end

end

error_msg='Max error: %0.4f';
fprintf(error_msg, error);

p3=plot(x,y2);
title("h=0.01: ");

hold on
p4=plot(x,f_array);
hold off

legend([p3 p4],{'Approximation','Function'});


%Re-do with h=0.0001: 
figure();
y2=[];
x=linspace(0,5,50000);

error=0;
for i=1:length(x)
    if x(i)==0
        y2(i)=1;
        f_array(i)=1;
    else
        y2(i)=y2(i-1) + 0.0001*((-2)*(y2(i-1)) + 2 - exp((-4)*x(i)));
        f_array(i)=f(x(i));
    end
    if abs(f_array(i)-y2(i))>error

        error=abs(f_array(i)-y2(i));

    end
end
error_msg='Max error: %0.4f';
fprintf(error_msg, error);

p4=plot(x,y2);
title("h=0.001: ");

hold on
p5=plot(x,f_array);
hold off

legend([p4 p5],{'Approximation','Function'});




g = @(t) 1/(1-t);

%Problem 2:
g_array=[];

t=linspace(0,0.99,99);


figure();
y4=[];

error=0;
for i=1:length(t)
    if t(i)==0
        y4(i)=1;
        g_array(i)=1;
    else
        y4(i)=y4(i-1) + 0.01*(y4(i-1))^2;
        g_array(i)=g(t(i));
    end
    if abs(f_array(i)-y(i))>error

        error=abs(g_array(i)-y(i));

    end
end
error_msg='Max error: %0.4f';
fprintf(error_msg, error);

p6=plot(t,y4);
title("g(t), h=0.01: ");

hold on
p7=plot(t,g_array);
hold off

legend([p6 p7],{'Approximation','Function'});

%Re-do with 0.001:

t=linspace(0,0.99,990);
g_array=[];


figure();
y4=[];
error=0;
for i=1:length(t)
    if t(i)==0
        y4(i)=1;
        g_array(i)=1;
    else
        y4(i)=y4(i-1) + 0.001*(y4(i-1)^2);
        g_array(i)=g(t(i));

    end
    if abs(f_array(i)-y4(i))>error

        error=abs(g_array(i)-y4(i));

    end
end
error_msg='Max error: %0.4f';
fprintf(error_msg, error);

p6=plot(t,y4);
title("g(t), h=0.001: ");

hold on
p7=plot(t,g_array);
hold off

legend([p6 p7],{'Approximation','Function'});


%Re-do with 0.0001:
g_array=[];

t=linspace(0,0.99,9900);


figure();
y4=[];

error=0;
for i=1:length(t)
    if t(i)==0
        y4(i)=1;
        g_array(i)=1;
    else
        y4(i)=y4(i-1) + 0.0001*y4(i-1)^2;
        g_array(i)=g(t(i));
    end
    if abs(f_array(i)-y4(i))>error

        error=abs(g_array(i)-y4(i));

    end
end
error_msg='Max error: %0.4f';
fprintf(error_msg, error);

p6=plot(t,y4);
title("g(t), h=0.0001: ");

hold on
p7=plot(t,g_array);
hold off

legend([p6 p7],{'Approximation','Function'});



lab_responses="\n I noticed that rather than having a consistent shrinking of the error difference, there was a spike in max error at h=0.001 and then a drastic shrinking as h=0.0001. This shows that it is possible to get a significantly low error margin with a slightly larger step size. Because smaller step size estimates require more processing power/effort, it is sometimes better to accept a larger error margin for less precision. ";
fprintf(lab_responses);


y(i)=y(i-1)+0.01*f(y(i-1)); 
