
syms F phidot u xcdot D gamma beta alpha phi

gamma = 2;%uncomment for 2e
alpha = 1;%uncomment for 2e
beta = 1;%uncomment for 2e
D = 1; %uncomment for 2e
mu = 3;%uncomment for 2e

xdd = [F-beta*(phidot^2)*sin(phi)-(mu*xcdot); D*sin(phi)];
pdd = [gamma -(beta*cos(phi)); -(beta*cos(phi)) alpha];
q = inv(pdd)*xdd;


q

