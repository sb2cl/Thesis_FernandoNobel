function [out] = simulate_host(p_base, t_inc)

t_inc = t_inc*60;

%% Intial condition simulation.

% Init model.
m = antithetic_controller();

% Parameters.
p = mergeStruct(p_base, m.p);
% p.cell__x1__omega_max = 0;
p.cell____burden = 0;
% p.cell__x1_openloop  = 0;

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [0 1e+9];

[~,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,m.x0,opt);

%% 1º part of simulation (No perturbation).

% Init model.
m = antithetic_controller();

% Parameters.
p = mergeStruct(p_base, m.p);
p.cell____burden = 0;

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [0 t_inc(1)];

[t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x(end,:),opt);
out = m.simout2struct(t,x,p);

%% 2º part of simulation (Perturbation).

% Init model.
m = antithetic_controller();

% Parameters.
p = mergeStruct(p_base, m.p);
p.cell____burden = 1;

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [t_inc(1) t_inc(1)+t_inc(2)];

[t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x(end,:),opt);
out = concatStruct(out, m.simout2struct(t,x,p));

%% 3º part of simulation (Recover from perturbation).

% Init model.
m = antithetic_controller();

% Parameters.
p = mergeStruct(p_base, m.p);
p.cell____burden = 0;

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [t_inc(1)+t_inc(2) t_inc(1)+t_inc(2)+t_inc(3)];

[t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x(end,:),opt);
out = concatStruct(out, m.simout2struct(t,x,p));

out.t = out.t/60;

end

