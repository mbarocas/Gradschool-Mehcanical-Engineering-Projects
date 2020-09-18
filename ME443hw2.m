clc; clear all; close all;
%% Universal Constants
kB= 1.38E-23; %J/k (Boltzman Constant)
hbar = 6.626E-34/(2*pi); %Js (Reduced Planck's Constant)
h = 6.626E-34; %Js (Planck's Constant)
c = 3E8; %m/s (Speed of Light)
sigma = 5.67E-8; %W m^-2 K^-4 (Stefan-Boltzman Constant - Photonic Radiation)
%% Question 1
% T = 500; %K
% Cphoton = (4*pi^2*kB^4*T^3)/(15*hbar^3*c^3);
% C = 16*sigma*T^3/c;
%% Question 2b
y = 0.3; %m^2/s
t  = 0;
x  = linspace(-15, 15, 1000);
k1 = 5; %m
k2 = 5.4; %m
kavg = 0.5*(k1+k2);
w1 = y*k1^2;
w2 = y*k2^2;
E1 = cos(w1*t-k1*x);
E2 = cos(w2*t-k2*x);
E  = E1+E2; %E Combined
figure
subplot(2,2,1);plot(x,E1,'r');xlabel('Distance (m)');ylabel('E-Field Wavepacket')
axis tight
title('E1')
subplot(2,2,2);plot(x,E2,'b');xlabel('Distance (m)');ylabel('E-Field Wavepacket')
axis tight
title('E2')
subplot(2,2,3);plot(x,E,'g');xlabel('Distance (m)');ylabel('E-Field Wavepacket')
%Xc Location
indexmin = find(min(E) == E); 
xmin = x(indexmin); 
Emin = E(indexmin);
indexmax = find(max(E) == E);
xmax = x(indexmax);
Emax = E(indexmax);
z = vline(xmax,'r','Xc = 0');
axis tight
title('E = E1+ E2')
subplot(2,2,4);hold on; plot(x,E1);plot(x,E2);plot(x,E,'g');xlabel('Distance (m)');ylabel('E-Field Wavepacket');legend('E1', 'E2', 'E')
axis tight
title('E1, E2, and E Superimposed')
%% Question 2c
% Now plot E(x) for three additional times, t = 1 s, t = 2 s, and t = 3 s. By looking at each plot,
% estimate xc, the location of the center of the wavepacket at time t. For example, at t = 0
% s, xc = 0 m. Since you are estimating xc by eye, do not worry about being perfectly accurate. 
y = 0.3; %m^2/s
t1 = 1;
t2 = 2;
t3 = 3;
x  = linspace(-5, 15, 1000);
k1 = 5; %m
k2 = 5.4; %m
w1 = y*k1^2;
w2 = y*k2^2;
E11 = cos(w1*t1-k1*x);
E12 = cos(w1*t2-k1*x);
E13 = cos(w1*t3-k1*x);
E21 = cos(w2*t1-k2*x);
E22 = cos(w2*t2-k2*x);
E23 = cos(w2*t3-k2*x);
Et1 = E11+E21;
Et2 = E12+E22;
Et3 = E13+E23;

figure
% E @ T=1
subplot(3,1,1);plot(x,Et1,'k');xlabel('Distance (m)');ylabel('E-Field Wavepacket')
%Xc Location
indexmin = find(min(Et1) == Et1); 
xmin = x(indexmin); 
Et1min = Et1(indexmin);
indexmax = find(max(Et1) == Et1);
xmax1 = x(indexmax);
Et1max = Et1(indexmax);
z1 = vline(xmax1,'r','Xc = 2.7678');
axis tight
title('E at t = 1')

% E @ T=2
subplot(3,1,2);plot(x,Et2,'k');xlabel('Distance (m)');ylabel('E-Field Wavepacket')
%Xc Location
indexmin = find(min(Et2) == Et2); 
xmin = x(indexmin); 
Et2min = Et2(indexmin);
indexmax = find(max(Et2) == Et2);
xmax2 = x(indexmax);
Et2max = Et2(indexmax);
z2 = vline(xmax2,'r','Xc = 6.7518');
axis tight
title('E at t = 2')

% E @ T=3
subplot(3,1,3);plot(x,Et3,'k');xlabel('Distance (m)');ylabel('E-Field Wavepacket')
indexmin = find(min(Et3) == Et3); 
xmin = x(indexmin); 
Et3min = Et3(indexmin);
indexmax = find(max(Et3) == Et3);
xmax3 = x(indexmax);
Et3max = Et3(indexmax);
z3 = vline(xmax3,'r','Xc = 9.5145');
axis tight
title('E at t = 3')

%% Question 2d

% Plot Xc as a function of T
tx = [0,1,2,3];
Xc = [0, xmax1, xmax2, xmax3];
figure
hold on
plot(tx,Xc,'k.');xlabel('Time');ylabel('Location of Center of Wavepacket');

% Comparison to Xc,Phase and Xc,Group
k1 = 5; %m
k2 = 5.4; %m
kavg = 0.5*(k1+k2);
w = y*(kavg)^2;
vphase = y*kavg;
vg = 2*y*kavg;
xcphase = vphase*tx;
xcg = vg*tx;
legend(Xc,xcphase,xcg)
plot(tx,xcphase,'r'); % Phase Velocity
plot(tx,xcg,'b'); % Group Velocity
%% Question 4
N=2;
Veff = 6000;
puc = 2.5E28;
kd = (6*pi^2*puc)^(1/3);
wd = Veff*kd;
thetad = wd*hbar/kB;
we = 8.8e13;
T = 200;
x = (hbar*we)/(kB*T)
8.28e-23*(x^2*exp(x))/((exp(x)-1)^2)