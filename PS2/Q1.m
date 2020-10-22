%% Q3
clear all;
close all;
%Ref HW1 Q2h
tspan = 0:0.01:200;

gamma = 2;
alpha = 1;
beta = 1;
D = 1;
mu = 3;

R = 10;
Q = [1000 0 0 0 ; 0 5 0 0 ; 0 0 1 0 ; 0 0 0 5];
x_00 = [0; 0; 0; 0];        % initial condition for 2g

A = [0,0,1,0;0,0,0,1;0,1,-3,0;0,2,-3,0]              
B = [0;0;1;1]
C = [39,0,0,0];
D = [0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0];

[K,S,e] = lqr(A,B,Q,R)
  
[tout,xout] = ode45(@(t,x) trackingController(t,x,K,C,A,B), tspan, x_00);

figure
%plot(tout,xout)
grid
title('Tracking Controller State vs Time')

%yout = C.*xout;

tf=200;
T = 0.01;
yd1= 20*ones(1,(tf/4)/T);
yd2= -20*ones(1,(tf/4)/T);
yd = cat(2, yd1,yd2,yd1,yd2);

figure
plot(yd)
hold on
%plot(yout)
grid
title('Desired vs Actual Output over time')
legend('Desired','Actual')
hold off

% OBSERVER PLACEMENT 
sys3 = ss(A,B,C,0)
ob = obsv(sys3)
observability = rank(ob)
sysPoles = eig(A) %system poles


A_transp = transpose(A);
C_transp = transpose(C);

%Kc from LQR :
Kc = [-9.99999999999901, 34.1720722633026, -28.5448907859479, 33.1569959170274]

%define observer poles as >5x dominant system poles
obsPole = [-20, -30, -40, -50]


L = place(A_transp, C_transp, obsPole)
Ko = transpose(L)


xhatdot = A*x + B*u - Ko(y - yhat)

%%%%%%%%%%%%%%%%


tf = 200;
step = 0.01;
numIteration = ceil (tf/step);
x = zeros(1,numIteration);
xhat = zeros(1,numIteration);
xdot = zeros(1,numIteration);
xhatdot = zeros(1,numIteration);
zdot = [A, -B*Kc; Ko*C, A-B*Kc-Ko*C]

A_fb = zdot;
disp('Q5.2 resultant Poles of A_fb:')
disp(eig(A_fb));


for k = 1:numIteration-1
    k1 = T*(A_fb*x(:,k));
    k2 = T*(A_fb*(x(:,k)+0.5*k1));
    k3 = T*(A_fb*(x(:,k)+0.5*k2));
    k4 = T*(A_fb*(x(:,k)+k3));
    x(:,k+1) = x(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
    t(k+1) = t(k) + T;
    z
end 
    

function xdot = trackingController(t,x,K,C,A,B);

x_2 = x(2);%phi
x_3 = x(3);%xcdot
x_4 = x(4);%phidot
yd=0;
if t<50
        yd = 20;
    end
    if t>=50 && t<100
        yd = -20;
    end
    if t>=100 && t<150
        yd = 20;
    end
    if t>=150 && t<200
        yd = -20;
    end

ABK= A-B*K;
ABKinv =inv(ABK);

v_partial = -(C*(ABKinv)*B)^(-1);
v = v_partial * yd;
F = v-K*x;

xcddot=(sin(x_2)*x_4^2 - F + 3*x_3)/(cos(x_2)^2 - 2) - (cos(x_2)*sin(x_2))/(cos(x_2)^2 - 2)
phiddot=(cos(x_2)*(sin(x_2)*x_4^2 - F + 3*x_3))/(cos(x_2)^2 - 2) - (2*sin(x_2))/(cos(x_2)^2 - 2)
xdot = [x_3;x_4;xcddot;phiddot];
end