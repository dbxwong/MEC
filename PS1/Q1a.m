
%Q1(a)(b)
A = [0,1,0;0,0,1;1,5,7]                         % Define System Matrices
B = [1;0;0]
C = [0,1,3]
x_0 = [0;1;0]
D = 0
ts = 2
eigA = eig(A)

Co = ctrb(A,B)
%unco = length(A) - rank(Co)

%Q1(c)
                                  
N = 50;                                    % Length Of Simulation (Samples)
t = linspace(0, 2, N);                     % Time Vector
u = ones(2, N);                            % System Input
for k1 = 1:length(t);
    y(:,k1) = C*expm(A*t(k1))*x_0.*u(:,k1);  % Evaluate System At Each Time & Input
end
%figure(1)
%plot(t, y)
%title('Output of unforced system for t=[0,2]')
%xlabel('Time (t)')
%ylabel('Output (y)')
%grid

%%%

%Q1(d)
p = [-1+i,-1-i,-2]
K = place(A,B,p)

%Q1(e)
Z_eigenvalues = eig(A,'matrix')
N = 50;                                    % Length Of Simulation (Samples)
t = linspace(0, 10, N);                     % Time Vector
u = ones(2, N);                            % System Input
for k1 = 1:length(t);
    z(:,k1) = C*expm(Z_eigenvalues*t(k1))*x_0.*u(:,k1);  % Evaluate System At Each Time & Input
end
%figure(2)
%plot(t, z)
%title('System output for t=[0,10]')
%xlabel('Time (t)')
%ylabel('Output (z)')
%grid


syms F phidot u xcdot D gamma beta alpha phi

gamma = 2;%uncomment for 2e
alpha = 1;%uncomment for 2e
beta = 1;%uncomment for 2e
D = 1; %uncomment for 2e
u = 3;%uncomment for 2e

xdd = [F-beta*(phidot^2)*sin(phi)-(u*xcdot); D*sin(phi)];
pdd = [gamma -(beta*cos(phi)); -(beta*cos(phi)) alpha];
q = inv(pdd)*xdd;
q




%Q2(c) this is for the linearized system
clear
A2 = [0,0,1,0;0,0,0,1;0,1,-3,0;0,2,-3,0]             
B2 = [0;0;1;1]
eigA2 = eig(A2)
Co2 = (ctrb(A2,B2));
unco2 = length(A2) - rank(Co2)


%2(d)
R = 10;
Q = [1,0,0,0 ; 0,5,0,0 ; 0,0,1,0 ; 0,0,0,5];
K2 = lqr(A2,B2,Q,R) %for linearized system
tspan = [0:0.01:30];

%initial conditions
x_01 = [0; 0.1; 0; 0];
x_02 = [0; 0.5; 0; 0];
x_03 = [0; 1.0886;0;0];
x_04 = [0;1.1;0;0];



%use anon function - A2,B2,K2 for linearized system
dx=@(t,x)((A2-B2*K2)*x);

%use ODE45 to solve the ode - x_01
[t,x]=ode45(dx,tspan,x_01);
plot(t,x)
grid
title('State of Linearized System under Feedback Control(ODE solved), IC:x_0_1')
xlabel('Time (t)')
ylabel('State (x(t))')

%use ODE45 to solve the ode - x_02
[t,x]=ode45(dx,tspan,x_02);
plot(t,x)
grid
title('State of Linearized System under Feedback Control(ODE solved), IC:x_0_2')
xlabel('Time (t)')
ylabel('State (x(t))')

%use ODE45 to solve the ode - x_03
[t,x]=ode45(dx,tspan,x_03);
plot(t,x)
grid
title('State of Linearized System under Feedback Control(ODE solved), IC:x_0_3')
xlabel('Time (t)')
ylabel('State (x(t))')

%use ODE45 to solve the ode - x_04
[t,x]=ode45(dx,tspan,x_04);
plot(t,x)
grid
title('State of Linearized System under Feedback Control(ODE solved), IC:x_0_4')
xlabel('Time (t)')
ylabel('State (x(t))')
%} 

%Q2(e)
close all
clear all
x_01 = [0; 0.1; 0; 0];
x_02 = [0; 0.5; 0; 0];
x_03 = [0; 1.0886;0;0];
x_04 = [0;1.1;0;0];

gamma = 2;
alpha = 1;
beta = 1;
D = 1;
mu = 3;

R = 10;
Q = [1 0 0 0 ; 0 5 0 0 ; 0 0 1 0 ; 0 0 0 5];

A2 = [0,0,1,0;0,0,0,1;0,1,-3,0;0,2,-3,0]              
B2 = [0;0;1;1]

K2 = lqr(A2,B2,Q,R) %for linearized system
tspan = [0:0.01:30]


[tout,xout] = ode45(@(t,x) nonlinear_sys(x,K2,alpha,beta,gamma,D,mu),tspan, x_01);
figure 
plot(tout,xout)
grid
title('State of Non-linearized System under Feedback Control(ODE solved), IC:x_0_1')
xlabel('Time (t)')
ylabel('State (x(t))')

[tout,xout] = ode45(@(t,x) nonlinear_sys(x,K2,alpha,beta,gamma,D,mu),tspan, x_02);
figure 
plot(tout,xout)
grid
title('State of Non-linearized System under Feedback Control(ODE solved), IC:x_0_2')
xlabel('Time (t)')
ylabel('State (x(t))')

[tout,xout] = ode45(@(t,x) nonlinear_sys(x,K2,alpha,beta,gamma,D,mu),tspan, x_03);
figure 
plot(tout,xout)
grid
title('State of Non-linearized System under Feedback Control(ODE solved), IC:x_0_3')
xlabel('Time (t)')
ylabel('State (x(t))')

[tout,xout] = ode45(@(t,x) nonlinear_sys(x,K2,alpha,beta,gamma,D,mu),tspan, x_04);
figure 
plot(tout,xout)
grid
title('State of Non-linearized System under Feedback Control(ODE solved), IC:x_0_4')
xlabel('Time (t)')
ylabel('State (x(t))')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%anonymous function for 2e non-linear case
function xdot = nonlinear_sys(x,K2,alpha,beta,gamma,D,mu)

F = -K2*x;
x_2 = x(2);%phi
x_3 = x(3);%xcdot
x_4 = x(4);%phidot

xcddot = (alpha*F+beta*D * cos(x_2)*sin(x_2) - alpha*beta*x_4^2*sin(x_2) - alpha*mu*x_3) / (alpha*gamma - beta^2*cos(x_2)^2);
phiddot = (xcddot*beta&cos(x_2) + D*sin(x_2))/alpha;

xdot = [x_3;x_4;xcddot;phiddot];
end

