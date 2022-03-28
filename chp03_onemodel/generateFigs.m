addpath('./build');

%% ex01_simple_gene_expression

clear all;
close all;

% Init model.
m = ex01_simple_gene_expression();

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [m.opts.t_init m.opts.t_end];

[t,x] = ode15s(@(t,x) m.ode(t,x,m.p),tspan,m.x0,opt);
out = m.simout2struct(t,x,m.p);

% Plot result.
fig = figure('units','centimeters','position',[0,0,6,6]);

hold on;
plot(out.t, out.mRNA, 'LineWidth',1.4);
plot(out.t, out.protein, 'LineWidth',1.4);

grid on;
xticks([0:2:10]);
ylabel('Concentration [a.u.]', 'interpreter', 'latex');
xlabel('Time [a.u.]', 'interpreter', 'latex');
legend('mRNA','protein','Location','southeast', 'interpreter', 'latex');

print(fig,'./figs/ex01_simple_gene_expression.eps','-depsc');

%% ex05_protein_induced

clear all;
close all;

% Init model.
m = ex05_protein_induced();

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [m.opts.t_init m.opts.t_end];

[t,x] = ode15s(@(t,x) m.ode(t,x,m.p),tspan,m.x0,opt);
out = m.simout2struct(t,x,m.p);

% Plot result.
fig = figure('units','centimeters','position',[0,0,6,6]);

hold on;
plot(out.t, out.A__mRNA, 'LineWidth',1.4);
plot(out.t, out.A__protein, 'LineWidth',1.4);
plot(out.t, out.B__mRNA, 'LineWidth',1.4);
plot(out.t, out.B__protein, 'LineWidth',1.4);

grid on;
xticks([0:2:10]);
ylabel('Concentration [a.u.]', 'interpreter', 'latex');
xlabel('Time [a.u.]', 'interpreter', 'latex');
legend('A.mRNA','A.protein','B.mRNA','B.protein','Location','southeast', 'interpreter', 'latex');

print(fig,'./figs/ex05_protein_induced.eps','-depsc');

%% ex06_antithetic_controller

clear all;
close all;

% Init model.
m = ex06_antithetic_controller();

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [m.opts.t_init 30];

[t,x] = ode15s(@(t,x) m.ode(t,x,m.p),tspan,m.x0,opt);
out = m.simout2struct(t,x,m.p);

% Plot result.
fig = figure('units','centimeters','position',[0,0,6,6]);

hold on;
plot(out.t, out.circuit__z1__protein, 'LineWidth',1.4);
plot(out.t, out.circuit__z2__protein, 'LineWidth',1.4);
plot(out.t, out.circuit__x__protein, 'LineWidth',1.4);

grid on;
% xticks([0:2:10]);
ylabel('Concentration [a.u.]', 'interpreter', 'latex');
xlabel('Time [a.u.]', 'interpreter', 'latex');
legend('z1.protein','z2.protein','x.protein','Location','best', 'interpreter', 'latex');

print(fig,'./figs/ex06_antithetic_controller.eps','-depsc');
