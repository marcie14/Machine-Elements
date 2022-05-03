% M E 338
% Calculations
% Group 2

%% Prepare Workspace
clear variables
close all
clc 
%% Global Variables:
% General Information
mass = 3; % Mass of car and components [kg]
Cd = 1; % Drag coefficient
g = 9.81; % Acceleration due to gravity [m/s^2]
pi = 3.14159; 
Fn = mass * g; % Normal Force [N]
rho = 1.2; % Air density [kg/m^3]

% Wheel & Axle Information
d = .00635; % diameter Axle [m] (1/4 inch)
I = pi * d^4 / 64; % Moment of Inertia Axles
A_axle = pi * d^2 / 4; % cross-section Axle
L_axle = 0.23; % length of the axle [m]
rwheels = 0.0635 / 2; % Radius of wheels [m]
distaxles = 0.2286; % Distance between axles [m]
distwheels = 0.168; % Distance between wheels on same axle [m]
Sya = 310; % Yield Strength Axle [MPa]
Ea = 2E9 / 10^6; % Elastic Modulus of 1045 Carbon Steel [MPa]
Crr = 0.02; % Rolling Resistance of car tire on asphalt
StaticF = 0.85; % static friction coefficient
KineticF = 0.85; % kinetic friction coefficient
a = .025; % distance between F wheels to F bearing

% Chassis Information
wood_width = .00635; % width of chassis [m]
wood_front_depth = .2413; % depth of chassis [m]
A_chassis = wood_width * wood_front_depth; % Area of Chassis [m^2]
L_Chassis = 0.305; % Length of chassis [m]
clearance = 0.025; % Ground clearance [m]
height = 0.05; % Height of front cross-section [m]
width = 0.2; % Width of front cross-section [m]
Sutc = 30; % UTS Chassis Baltic Birch Plywood [MPa]
Succ = 36; % UCS Chassis Baltic Birch Plywood [MPa]
Ec = 9E9 / 10^3; % Elastic Modulus of Baltic Birch Plywood [MPa]

% Gear Drive Information
TR = 0.16; % Transmission Ratio
Efficiency = 0.8; % Gear Train Efficiency
rpm_motor = 22500; % [rpm]

% Impact Forces 
nu = 1; % Correction factor due to energy dissipation 
Test_N = 100; % Solidworks test load [N]

disp("-------------------------------------------------------------------------");
%% Max Velocity and Acceleration
% Velocity
rpm_wheel = rpm_motor * TR; % RPM of wheel given RPM motor and TR [rpm]
V_wheel = 2 * pi * rwheels * rpm_wheel / 60; % Max Velocity [m/s]

% Forces and Torques
Fr = Crr * Fn + .5 * A_chassis * rho * V_wheel^2; % Resistance Force = Rolling + Drag
T_wheelsV = Fr * rwheels; % Torque at Wheels [mN*m] 
T_motor = T_wheelsV / Efficiency; % Torque at motor [mN*m] 

% Acceleration
T_wheelsA = 75 * T_motor / 12; % Torque at Wheels [N*m]
MaxF_wheels = T_wheelsA / rwheels; % Max Force at Wheels [N]
MaxF_noslip = Fn * (StaticF + Crr); % Max Force at Wheels (no slip) [N]
MaxA_noslip = MaxF_noslip / mass; % Max Acceleration at Wheels (no slip) [m/s^2]

disp("Maximum Velocity: " + V_wheel + " m/s");
disp("Maximum Force: " + MaxF_wheels + " N");
disp("Maximum Force (no slipping): " + MaxF_noslip + " N");
disp("Maximum Acceleration (no slipping): " + MaxA_noslip + " m/s^2");
disp("-------------------------------------------------------------------------");

%% Reactions
disp("Chassis and Axle Reaction Forces");
% Reactions at Max Velocity
A_chassis = width * height; % Front cross-sectional area [m^2]
Fd = (1/2) * Cd * rho * V_wheel^2 * A_chassis; % Drag force [N]
Rf_Vm = mass * g / 2 - Fd * (0.5 * height + clearance) / distaxles; % Vertical reaction at front wheels [N]
Rr_Vm = mass * g / 2 + Fd * (0.5 * height + clearance) / distaxles; % Vertical reaction at rear wheels [N]
disp("Max Velocity");
disp("     Vertical Reaction at Front Wheels: " + Rf_Vm + " N");
disp("     Vertical Reaction at Rear Wheels: " + Rr_Vm + " N");

% Reactions at Max Acceleration
Ft = mass * MaxA_noslip; % Force of friction accelerating car [N]
Rf_Am = mass * g / 2 - Ft * rwheels / distaxles; % Vertical reaction at front wheels [N]
Rr_Am = mass * g /2 + Ft * rwheels / distaxles; % Vertical reaction at rear wheels [N]
disp("Max Acceleration");
disp("     Vertical Reaction at Front Wheels: " + Rf_Am + " N");
disp("     Vertical Reaction at Rear Wheels: " + Rr_Am + " N");
disp("-------------------------------------------------------------------------");

%% Impact Forces
% Impact Force
disp("Chassis and Axle Impact Forces");
Numer = nu * mass * Ec * A_chassis * 48 * Ea * I;
Denom = 4 * Ec * A_chassis * a^3 - 3 * Ec * A_chassis * a * distaxles^2 + 48 * Ea * L_Chassis * I ;
Fi = (V_wheel * sqrt(-1 *Numer / Denom)); % Impact Force [N]
disp("     Impact Force: " + Fi + " N");

% Stiffness
deflection_axle = Test_N * A_axle * (4 * A_axle^2 - 3 * L_axle^2) / (24 * Ea * I);
deflection_chassis = 2.457E-5; % Solidworks deflection with a 100N load [m]
k_axle = Test_N / deflection_axle; % stiffness of axle
k_chassis = Test_N / deflection_chassis; % stiffness of chassis
k_total = ((k_axle * k_chassis) / (k_axle + k_chassis)); % total stiffness
disp("     Stiffness: " + k_total + " kN/m");

disp("-------------------------------------------------------------------------");

%% Impact Stresses
disp("Impact Stresses");

% Axle
disp("Axle");
MaxMomA = Rr_Am / 2 * a; % Max Moment Axles [N*m]
Va = Rr_Am / 2; % Shear [N]
sigxa = MaxMomA / (I * 10^6); % Sigmax Axle ** Was hard coded to 300 ** [MPa]
sigya = 0 / 10^6; % Sigmay Axle [MPa]
txya = 4 * Va / (3 * A_axle * 10^6); % Max shear Axles [MPa]
disp("     Impact Stress (Sigma x): " + sigxa + " MPa");
disp("     Impact Shear Stress (Tau xy): " + txya + " MPa");

% Chassis
disp("Chassis");
sigAc = 0; % SigmaA Chassis [MPa]
sigBc = -Fi / (A_chassis * 10^6); % SigmaB Chassis ** -3.53 ** [MPa]
txyc = abs((sigAc - sigBc) / 2 ); % Tauxy Chassis [MPa] ** 1.765**
disp("     Impact Stress (Sigma A): " + sigAc + " MPa");
disp("     Impact Stress (Sigma B): " + sigBc + " MPa");
disp("     Impact Shear Stress (Tau 1,2): " + txyc + " MPa");
disp("-------------------------------------------------------------------------");

%% Axle Safety Factor
disp("Axle Safety Factors");

% Distortion Energy Theory (von Mises)
sig1a = (sigxa + sigya)/2 + sqrt(((sigxa - sigya)/2)^2 + txya^2);
sig3a = (sigxa + sigya)/2 - sqrt(((sigxa-sigya)/2)^2 + txya^2);

% if either sigma is 0:
if (sigxa == 0 || sigya == 0)
    sigprime = sqrt(sigxa^2 - sigxa*sigya + sigya^2 + 3*txya^2);
    %disp("**Alternate sigma prime equation used since one sigma is 0**");
else % otherwise:
    sigprime = sqrt(sig1a^2 - sig1a*sig3a + sig3a^2);
end

n1 = Sya/sigprime;

disp("     Axle Distortion Energy Safety Factor: " + n1);

% Maximum Shear Stress Theory
  
% if both principal stresses are positive, shear is maximized when third at 0 is included
if (sig1a > 0 && sig3a > 0) 
    sig3a = 0;
    %disp("**sig3 = 0 for maximizing shear stress since both principal stresses are positive**");
end

tmax = abs((sig1a - sig3a))/2;
n2 = (Sya/2)/tmax;

disp("     Axle Maximum Shear Stress Safety Factor (More Conservative): " + n2);

%% Chassis Safety Factors
disp("Chassis Safety Factors");

% Coulumb-Mohr & Modified Mohr
if ((sigAc >= sigBc) && (sigBc >= 0))
    n3 = Sutc/sigAc; % Coulomb-Mohr
    n4 = n3; % Modified-Mohr
elseif ((sigAc >= 0) && (0 >= sigBc))
    n3 = (Sutc*Succ) / (Succ*sigAc - Sutc*sigBc); % Coulomb-Mohr
    if ((abs(sigBc / sigAc) <= 1) && sigAc ~= 0)
        n4 = Sutc / sigAc; % Modified-Mohr
    else
        n4 = (Sutc*Succ) / (Succ*sigAc + Sutc*(sigAc - sigBc)); % Modified-Mohr
    end
else
    n3 = Succ / sigBc; % Coulomb-Mohr
    n4 = n3; % Modified-Mohr   
end

disp("     Chassis Coulomb-Mohr Safety Factor: " + n3);
disp("     Chassis Modified-Mohr Safety Factor: " + n4);

disp("-------------------------------------------------------------------------");

