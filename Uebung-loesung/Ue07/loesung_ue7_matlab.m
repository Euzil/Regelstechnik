% Übung 7 - 04.06.2021

clear variables
close all
clc

% Vorbereitung
s = tf('s'); % Matlab mitteilen, dass s Variable von TF/ÜF ist

%% A1 - Stabilitätsuntersuchung des geschlossenen Kreises (Nyquist)

% ÜF aufstellen
Ga = 1/((s+1)*(s+2));
Gb = s/((s+1)*(s+2));
Gc = 1/(s*(s+1)*(s+2));
Gd = 1/((s-1)*(s+2));
Ge = 1/((s+1)*(s-2));
Gf = 1/(s*(s-1)*(s+2));
Gg = (s+2)/(s*(s-1));

% Nyquistdiagramme erstellen lassen 
figure('Name','1a)')
nyquist(Ga)
figure('Name','1b)')
nyquist(Gb)
figure('Name','1c)')
nyquist(Gc)
figure('Name','1d)')
nyquist(Gd)
figure('Name','1e)')
nyquist(Ge)
figure('Name','1f)')
nyquist(Gf)
figure('Name','1g)')
nyquist(Gg)

%% A2 - Frequenzgang wichtiger Übertragungsfunktionen (Bode) 

% Variablen, um Formeln der Regler aufzustellen
% -> Variiert diese, um die unterschiedlichen Einflüsse auf die Magnitude
% und Phase zu sehen
kp = 15; % P-Anteil
Ti = 2; % I-Anteil
Td = 2; % D-Anteil
T1 = 10; % T für Lead
T2 = 1/10; % T für Lag

% Regler aufstellen
P = kp + 0*s;
PI = kp * (1 + Ti*1/s);
PD = kp * (1 + Td*s);
PID = kp * (1 + Td*s + Ti*1/s);
LEAD = (T1*s+1)/(0.5*T1*s+1);
LAG = (T2*s+1)/(1.5*T2*s+1);

% Bode Plots in Subfigures erstellen
figure('Name','Aufgabe 2)')
subplot(2,3,1)
bode(P)
title('P')
grid on
subplot(2,3,2)
bode(PI)
title('PI')
grid on
subplot(2,3,3)
bode(PD)
title('PD')
grid on
subplot(2,3,4)
bode(PID)
title('PID')
grid on
subplot(2,3,5)
bode(LEAD)
title('LEAD')
grid on
subplot(2,3,6)
bode(LAG)
title('LAG')
grid on

%% A3 Reglreentwurf im Frequenzbereich

% System Performance untersuchen
fig1 = figure('Name','Aufgabe 3 - Performance Übersicht');
hold on
% Durch zweiten Parameter wird ein Einheitssprung auf das unbekannte System
% gegeben (1 -> Einheitsprung). Beim ersten Parameter muss der Regler
% eingefügt werden (1 -> Kein Regler)
G_heim(1,1) 
legend('G')
% - Sehr langsamer Übergang (ca. 3s bis zum Steady State)
% - Große statische Regelabweichung (ca. 1/3 des Zielwerts)

% Bode Plot für das unbekannte System betrachten (nur ein Parameter -> 
% hier würde auch wieder der Regler eingefügt werden)
fig2 = figure('Name','Aufgabe 3 - Bode Stabilität Übersicht');
hold on
G_heim(1)
legend('G')

%% 1.
% Steady state gain (bode plot at s=0) of G_ss = -6dB = 1/2
% With a maximum steady state error of 0.05 it has to hold that:
% [1/(1+KG) < 0.05]	<=>   [K_ss*G_ss > 19]    <=>  [K_ss > 38]

K = 39;

figure(fig1)
G_heim(K,1)
legend('G','K*G')
figure(fig2)
G_heim(K)
legend('G','K*G')

%% 2. 
% Phasemargin >= 60°:
% Search for crossover frequency and lookup the current phasemargin (~30dB)
% Select controller:
%   a) PD-Compensator
%   b) LEAD-Compensator
%
% With LEAD: 
% - Select alpha. To do so, take a look at the table presented in lecture 8b.
%       +30 dB  --> alpha > 1/3 (Choose alpha a little smaller to compensate
%       the shift of the crossover to a higher frequency)
% - Select w_max and compute T=1/(w_max*sqrt(alpha)). (Experiment with
%   choosing w_max)

alpha = 1/5;
w_max = 12;
T = 1/(w_max*sqrt(alpha))
D_lead = K*((T*s+1)/(alpha*T*s+1));

figure(fig1)
G_heim(D_lead,1)
legend('G','K*G','K*D*G')
figure(fig2)
G_heim(D_lead)
legend('G','K*G','K*D*G')

%% 3.
% Exermine the effect of using an integrator instead of the static gain 
% to reduce the steady state error of (1.) any further.
%
% Use the very same controller as in (2.) but replace K=38 with K=1/s:

D = 1/s * 1/K * D_lead;

figure(fig1)
G_heim(D,1)
legend('G','K*G','K*D*G','1/s*D*G')
figure(fig2)
G_heim(D)
legend('G','K*G','K*D*G','1/s*D*G')

% ==> The controller with the integrator has a much smaller bandwith 
% (w_B < w_c)

%% A4
% Example for DC-engine (lecture)
%G = 1/(s*(s+1));

% more complexe system
G = (s+4)/(s*(s+1)*(s+10));

% with |G_ss| = 1, R = 1/s^2 and a steady state error of max 2%, we have 
% R/(1+KG) < 0.02	<=>     K > 50
K = 60;

% ==> Implement Lag-Compensation
% Using alpha = K
alpha=K; 
% Using w<7 rad/s as cut-off frequency for the zero of the lag compensator
T=10;
D = alpha* (T*s+1)/(alpha*T*s+1);

% plot results
figure('Name','Aufgabe 4 - Bode Stabilität Übersicht');

bode(1*G);
hold on
grid on
bode(K*G)
bode(D*G)
legend('G','K*G','D*G')