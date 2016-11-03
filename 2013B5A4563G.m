clear all %clears all variables
close all %closes all previously opened figures
clc
global K % Defined globally because of using in the function f1.m 
global Q  % Defined globally because of using in the function f1.m 
a = 0; %Starting value of Zbar
b = 1; %Last value of Zbar
N = 100; %No of Iterations
h = (b-a)/N;

K = 1; 
v = [0,-1,-5,-10] %vector defining different K
tol = 10^(-7); %tolerance


for k = 1:length(v)
    Q = v(k);
err = 4; %Some random value for error 
guess1 = 1.1; %Guess 1 for Secant Method
guess2 = 1.2; %Guess 2 for Secant Method
l = 0; % Boundary condition for lambda
l1 = 1; % Boundary condition for lambda'
l2 = guess1;  % Boundary condition for lambda'' this value is bought to zero by succesive iterations
for i = 1:N
    z(i) = (i-1)*h;
end
for j = 1:N
        lm = l  + l1*h/2; %calculating the value of lambda''(1) as a second guess for secant method
        l1m = l1 + l2*h/2; 
        l2m = l2 + f1(l,l1,l2)*h/2;
        
        lf = l + l1m*h; %calculating the value of lambda''(1) as a second guess for secant method (final timestep)
        l1f = l1 + l2m*h;
        l2f = l2 + f1(lm,l1m,l2m)*h;
        l=lf;
        l1=l1f;
        l2=l2f;
end
ans1=l2;
l1save=zeros(N,1);  %to save the values of lambda'
l0save=zeros(N,1); %to save the values of lambda
while err > tol
    l=0; %Boundary Conditions
    l1=1; %Boundary Conditions
    l2=guess2;
    for i = 1:N
        l1save(i)=l1;
        l0save(i)=l;
        lm = l  + l1*h/2; %First set of steps of RK-2
        l1m = l1 + l2*h/2;
        l2m = l2 + f1(l,l1,l2)*h/2;
        
        lf = l  + l1m*h;  %Second set of steps for RK-2
        l1f = l1 + l2m*h;
        l2f = l2 + f1(lm,l1m,l2m)*h;
        l=lf;
        l1=l1f;
        l2=l2f;
        
        end
    ans2=l2;
    err = abs(l2f); %calculation of error
    guess3 =  guess2 - ans2*(guess2-guess1)/(ans2-ans1); %Secant Method formula
    guess1=guess2; %Reassigning values
    guess2=guess3;
    ans1=ans2;
end
e_norm = 0.5*sqrt(2*((1-l1save).^2)+(1-1./l1save).^(2)); %Calculation of the norm
theta=-Q*l0save+1;  %Calcutating the non-dimesional temperature Eq 3.34 in the assignment
    figure(1)
    plot(z,l1save);
    xlabel('$\bar{Z}$','interpreter','Latex');
    ylabel('$\bar{\lambda\prime}$','interpreter','Latex');
    
    figure(2)
    plot(z,e_norm) %the values in the vector v must be changed
    xlabel('$\bar{Z}$','interpreter','Latex');
    ylabel('$ e(0,0,\bar(Z)) $' ,'interpreter','Latex');
    hold on
    grid on
    
   figure(3) 
   plot(z,theta);
   xlabel('$\bar{Z}$','interpreter','Latex'); %the values in the vector v must be changed
   ylabel('$\bar{\theta}$','interpreter','Latex');
   hold on
   grid on

hold on
end
l2
