close all; 
clear all;
clc; 

Equationsdata = fopen('SaturnV.dat','w');
dt = 1;
while dt>=1
% Earth and Launch Site Data
g0 = 9.81; % Initial gravity m/s^2
Re = 6378000; % meters
hs = 6000; % scale height
%h = altitude;
%rho = Rho; %kg/s^2

%LV Configuration 
CD = 0.06; % Coefficient of Drag
m1 = 131000; % mass kg
m2 = 36000;
m3 = 11000;
Sref1 = 451.229; % Planform area
Sref2 = 242.187;
Sref3 = 173.2263;
T1 = 3.3e+07; % Thrust 
T2 = 4.9e+06;
T3 = 1e+06;
m_dot1 = 250; % fuel flow rate
m_dot2 = 420;
m_dot3 = 420;
ff1 = T1/(g0*m_dot1);% mass flow rate
ff2 = T2/(g0*m_dot2);  
ff3 = T3/(g0*m_dot3);
L1 = 106.86; % Length meters
L2 = 63.39;
L3 = 44.11;
ve1 = 2.58e3; % exhaust velocity
ve2 = 4.13e3;
ve3 = 4.13e3; 

%Initial parameters 
m0 = 2951000;
pitch0 = 90;
pitch_rate = 0;
v0 = 0;
x0 = 0;
x_dot = 0;
h_dot = 0;
h0 = 0;
D = 0;
t0 = 0;
t = t0;
rho0 = 1.225;
a = T1/m0-g0; % only for the first timestep
fprintf(Equationsdata,'%f %f %f %f %f %f\n',t0,h0,pitch0,v0,a,m0);
t = t+dt;
    while t<153 % Stage 1 ignition
    m = m0 - ff1*t;
    v = v0 + a*t;
    h_dot = v*sind(pitch0);
    h = h0 + h_dot*t;
    rho = rho0*exp(-h/hs);
    D = CD*Sref1*rho.*0.5*v^2;
    g_local = g0/(1+h/Re).^2;
    Pitch_rate = -(g_local/v)+(v/Re+h0)*cosd(pitch0);
    pitch = pitch0 + Pitch_rate*t;
    x_dot = v*cosd(pitch0)*(Re/(Re+h));
    x = x0 + x_dot*t;
    a = ((T1-D)/m)-g_local*sind(pitch);
    fprintf(Equationsdata,'%f %f %f %f %f %f\n',t,h,pitch,v,a,m);
    t = t+dt;
    end

m = m0-m1; %Jettison stage 1

    while t<168  % no thrust 
    m = m;
    v = v + a*t;
    h_dot = v*sind(pitch0);
    h = h + h_dot*t;
    rho = rho*exp(-h/hs);
    D = CD*Sref1*rho.*0.5*v^2;
    g_local = g0/(1+h/Re).^2;
    Pitch_rate = -(g_local/v)+(v/Re+h)*cosd(pitch);
    pitch = pitch + Pitch_rate*t;
    x_dot = v*cosd(pitch)*(Re/(Re+h));
    x = x + x_dot*t;
    a = (-D/m)-g_local*sind(pitch);
    fprintf(Equationsdata,'%f %f %f %f %f %f\n',t,h,pitch,v,a,m);
    t = t+dt;
    %m_new = m0 + m_dot*dt;
    end

    while t<534.6 % ignite stage 2 
    m = m - ff2*t;
    v = v + a*t;
    h_dot = v*sind(pitch);
    h = h + h_dot*t;
    rho = rho*exp(-h/hs);
    D = CD*Sref1*rho.*0.5*v^2;
    g_local = g0/(1+h/Re).^2;
    Pitch_rate = -(g_local/v)+(v/Re+h)*cosd(pitch);
    pitch = pitch + Pitch_rate*t;
    x_dot = v*cosd(pitch)*(Re/(Re+h));
    x = x + x_dot*t;
    a = ((T2-D)/m)-g_local*sind(pitch);
    fprintf(Equationsdata,'%f %f %f %f %f %f\n',t,h,pitch,v,a,m);
    t = t+dt;
    end

m = m0-m1-m2; % jettison stage 2

    while t<537.8 % no thrust
    m = m;
    v = v + a*t;
    h_dot = v*sind(pitch);
    h = h + h_dot*t;
    rho = rho*exp(-h/hs);
    D = CD*Sref1*rho.*0.5*v^2;
    g_local = g0/(1+h/Re).^2;
    Pitch_rate = -(g_local/v)+(v/Re+h)*cosd(pitch);
    pitch = pitch + Pitch_rate*t;
    x_dot = v*cosd(pitch)*(Re/(Re+h));
    x = x + x_dot*t;
    a = (-D/m)-g_local*sind(pitch);
    fprintf(Equationsdata,'%f %f %f %f %f %f\n',t,h,pitch,v,a,m);
    t = t+dt;
    end

    while t<983.5 %ignite stage 3
    m = m - ff3*t;
    v = v + a*t;
    h_dot = v*sind(pitch);
    h = h + h_dot*t;
    rho = rho*exp(-h/hs);
    D = CD*Sref1*rho.*0.5*v^2;
    g_local = g0/(1+h/Re).^2;
    Pitch_rate = -(g_local/v)+(v/Re+h)*cosd(pitch);
    pitch = pitch + Pitch_rate*t;
    x_dot = v*cosd(pitch)*(Re/(Re+h));
    x = x + x_dot*t;
    a = ((T3-D)/m)-g_local*sind(pitch);
    fprintf(Equationsdata,'%f %f %f %f %f %f\n',t,h,pitch,v,a,m);
    t = t+dt;
    %m_new = m0 + m_dot*dt;
    end

m = m0-m1-m2-m3; %jettison stage 3

    while t<1000 % no thrust
    m = m;
    v = v + a*t;
    h_dot = v*sind(pitch0);
    h = h + h_dot*t;
    rho = rho*exp(-h/hs);
    D = CD*Sref1*rho.*0.5*v^2;
    g_local = g0/(1+h/Re).^2;
    Pitch_rate = -(g_local/v)+(v/Re+h)*cosd(pitch);
    pitch = pitch + Pitch_rate*t;
    x_dot = v*cosd(pitch)*(Re/(Re+h));
    x = x + x_dot*t;
    a = (-D/m)-g_local*sind(pitch);
    fprintf(Equationsdata,'%f %f %f %f %f %f\n',t,h,pitch,v,a,m);
    t = t+dt;
    end

    dt = dt*.1; 
end
fclose(Equationsdata);
load SaturnV.dat
% plots
figure
plot(SaturnV(:,1),SaturnV(:,4));
title('Acceleration')
ylabel('Acceleration [m/s^2]')
xlabel('time [s]')
hold on

figure
plot(SaturnV(:,1),SaturnV(:,2));
title('altitude')
ylabel('height [m]')
xlabel('time [s]')
hold on