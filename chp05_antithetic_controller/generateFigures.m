%% Base parameters for the simulation.

clear all;

% Base parameters for all the figures.
p_base = [];

K = 3; %4
p_base.cell__x1__omega_max = 5*K;
p_base.cell__x2__omega_max = 1*K;
p_base.cell__A__omega_max = 20;
p_base.cell__B__omega_max = 100;

p_base.cell__A__h = 2.5; % 5
p_base.cell____dilution = 1;
p_base.cell__gamma = 10;

% Time increment between simulation parts.
t_inc = [5 15 20];

%% Figures.

close all;

fig01(p_base, t_inc);
fig02(p_base, t_inc);
fig03(p_base, t_inc);
fig04(p_base, t_inc);
