%==========================================================================
%
% Advanced Topics in Control
%   Exercise - iPendel
%
% Institute for Electrical Engineering in Medicine
%   Dr.-Ing. Christian Hoffmann
%   Georg Männel
%
%   07/06/2017
%
%==========================================================================
%% Initialization
addpath(genpath('./'))
s  = tf('s');
Ts = 0.01;
%% Physical Plant Definition (linearized around the upright position)

% From Pendulum.c of the firmware:
% #define MB			0.138f		///< Pendulum body mass [kg]
% #define MW			0.094f		///< Pendulum wheel mass [kg]
% #define LB			0.0430f		///< Distance between ground and body center of mass [m]
% #define LW			0.0406f		///< Distance between ground and wheel center of mass [m]
% #define IB			9.85e-5f	///< Pendulum body inertia [kg.m^2]
% #define IW			4.20e-5f	///< Pendulum wheel inertia [kg.m^2]
% #define RA			0.015f		///< Distance between ground and accelerometer device [m]
% #define TRQ_B		    0.079f		///< Braking torque [N.m]
% #define KM			25.1e-3f	///< Motor torque constant [N.m/A]

u_max = 255e-3; % Maximum motor torque [N.m]

g        = 9.81;
mb       = 0.1380;
mw       = 0.0940;
lb       = 0.0430*1.5;
lw       = 0.0406*1.5;
Ib       = 9.85e-5;
Iw       = 4.20e-5;
ra       = 0.015;
trq_b    = 0.079;
km       = 25.1e-3;

% Total mass
M = mb*lb + mw*lw;
% Total inertia
I = Ib + mb*lb^2 + mw*lw^2;

% ddtheta = 1/I * (M*g*theta - u_m + u_f)
% ddphi   = 1/(I*Iw) * ((I + Iw)*(u_m - u_f) - 1/I * M*g*theta

 theta0 =  pi/4; %% Initial pendulum angle [rad]
dtheta0 = -pi/4; %% Initial pendulum angular velocity [rad/s]
  dphi0 =  1000; %% Initial motor speed [rad/s]
r_theta   = 0;     %% Go to upright equilibrium
r_dtheta  = 0;     %% Go to upright equilibrium with constant velocity
r_dphi    = dphi0; %% Keep initial motor speed

% Total mass
M = mb*lb + mw*lw;
% Total inertia
I = Ib + mb*lb^2 + mw*lw^2;

% Model linearized about the upright position
% x = [ theta, dtheta, dphi ];
A   = [   0          ,     1,      0 ; ...
          1/I * (M*g),     0,      0 ; ...
        - 1/I * (M*g),     0,      0 ];

% Gain of motor torque input
B  = [   0                  ;...
       - 1/I                ;...
         1/(I*Iw)*(I + Iw)  ];
    

C  = [1, 0, 0;...
      0, 0, 1];

D  = 0;

G   = ss(A, B , C, D );
ne = size(G, 1);
ny = size(G, 1);