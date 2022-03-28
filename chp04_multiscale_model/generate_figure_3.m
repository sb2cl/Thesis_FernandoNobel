clear all;
close all;

%% Step 1: Generate a 3-D array of values ku, kb, and N*omega (Gray dots).

% Limit values of gene expression space.
lim_ku    = [6 135];
lim_kb    = [3 15];
lim_omega = [100 350];

% Generate array with the values choosen to form the gene expression space.
ku    = arrayMinMaxN(lim_ku, 4);
kb    = arrayMinMaxN(lim_kb, 4);
% omega = [0.1 350];
% s = [0.1 3.6];

omega = [0.1 1 5 10 50 100 175 200 250 300 350];
s  = [0.1 0.2 0.4 0.6 0.8 1 1.5 2 2.5 3 3.6];

% Gene expression space.
input_gray = {};

% Fill the gene expression space.
for i = 1:length(ku)
    for j = 1:length(kb)
        for k = 1:length(omega)
            for l = 1:length(s)
                input_gray{i,j,k,l}.p.cell__p_A__k_u = ku(i);
                input_gray{i,j,k,l}.p.cell__p_A__k_b = kb(j);
                input_gray{i,j,k,l}.p.cell__p_A__omega = omega(k);
                input_gray{i,j,k,l}.p.bio__s = s(l);
            end
        end
    end
end

%% Step 2: Generate a small set of points varing N*omega and ku or kb (Color dots.)

% Generate array with the values choosen to form the gene expression space.
ku    = [6];
kb    = [3 9 15];
omega = [1 100 350];

% Gene expression space.
input_color = {};

% Fill the gene expression space.
for i = 1:length(ku)
    for j = 1:length(kb)
        for k = 1:length(omega)
            for l = 1:length(s)
                input_color{i,j,k,l}.p.cell__p_A__k_u = ku(i);
                input_color{i,j,k,l}.p.cell__p_A__k_b = kb(j);
                input_color{i,j,k,l}.p.cell__p_A__omega = omega(k);
                input_color{i,j,k,l}.p.bio__s = s(l);
            end
        end
    end
end

%% Step 3: Caclulate intial conditions for both sets with model_initialize.

% Gray.

% Init model.
m = model_initialize();

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [m.opts.t_init m.opts.t_end];

for i = 1:size(input_gray,1)
    for j = 1:size(input_gray,2)
        for k = 1:size(input_gray,3)
            for l = 1:size(input_gray,4)
                % Default simulation parameters.
                p = m.p;
                
                % Change the parameters of the gene expression space.
                p.cell__p_A__k_u   = input_gray{i,j,k,l}.p.cell__p_A__k_u;
                p.cell__p_A__k_b   = input_gray{i,j,k,l}.p.cell__p_A__k_b;
                p.cell__p_A__omega = input_gray{i,j,k,l}.p.cell__p_A__omega;
                p.s                = input_gray{i,j,k,l}.p.bio__s;
                
                % Simulate.
                [t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,m.x0,opt);
                out = m.simout2struct(t,x,p);
                
                % Create a struct with the steady-state info.
                x0 = [];
                x0.p_r__m  = out.cell__p_r__m(end);
                x0.p_nr__m = out.cell__p_nr__m(end);
                x0.mu      = out.cell__mu(end);
                x0.r       = out.cell__r(end);
                x0.p_A__m  = out.cell__p_A__m(end);
                
                % Save the steady-state as the intial condition for the next
                % simulations.
                input_gray{i,j,k,l}.x0 = x0;
            end
        end
    end
end

% Color.

% Init model.
m = model_initialize();

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [m.opts.t_init m.opts.t_end];

for i = 1:size(input_color,1)
    for j = 1:size(input_color,2)
        for k = 1:size(input_color,3)
            for l = 1:size(input_color,4)
            % Default simulation parameters.
            p = m.p;
            
            % Change the parameters of the gene expression space.
            p.cell__p_A__k_u   = input_color{i,j,k,l}.p.cell__p_A__k_u;
            p.cell__p_A__k_b   = input_color{i,j,k,l}.p.cell__p_A__k_b;
            p.cell__p_A__omega = input_color{i,j,k,l}.p.cell__p_A__omega;
            p.s                = input_color{i,j,k,l}.p.bio__s;
    
            % Simulate.
            [t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,m.x0,opt);
            out = m.simout2struct(t,x,p);
            
            % Create a struct with the steady-state info.
            x0.p_r__m  = out.cell__p_r__m(end);
            x0.p_nr__m = out.cell__p_nr__m(end);
            x0.mu      = out.cell__mu(end);
            x0.r       = out.cell__r(end);
            x0.p_A__m  = out.cell__p_A__m(end);
            
            % Save the steady-state as the intial condition for the next 
            % simulations.
            input_color{i,j,k,l}.x0 = x0;
            end
        end
    end
end

%% Step 4: Simulate both sets with model_fedbach.

% Gray.

output_gray = {};

m = model_fedbatch();

for i = 1:size(input_gray,1)
    for j = 1:size(input_gray,2)
        for k = 1:size(input_gray,3)
            for l = 1:size(input_gray,4)
                % Default simulation parameters.
                p = m.p;
                
                % Change the parameters of the gene expression space.
                p.cell__p_A__k_u   = input_gray{i,j,k,l}.p.cell__p_A__k_u;
                p.cell__p_A__k_b   = input_gray{i,j,k,l}.p.cell__p_A__k_b;
                p.cell__p_A__omega = input_gray{i,j,k,l}.p.cell__p_A__omega;
                
                % Default initial conditions.
                x0 = m.x0;
                
                % Use initial conditions of previous simulation.
                x0(1)  = input_gray{i,j,k,l}.x0.p_r__m;
                x0(2)  = input_gray{i,j,k,l}.x0.p_nr__m;
                x0(3)  = input_gray{i,j,k,l}.x0.mu;
                x0(4)  = input_gray{i,j,k,l}.x0.r;
                x0(5)  = input_gray{i,j,k,l}.x0.p_A__m;
                x0(10) = input_gray{i,j,k,l}.p.bio__s;
                
                % Solver options.
                opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
                opt = odeset(opt,'Mass',m.M);
                opt = odeset(opt,'Events',@(t,y) eventSubstrateDepletion(t,y,p,m));
                
                % Simulation time span.
                tspan = [0 1e9];
                
                % Simulate.
                [t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x0,opt);
                out = m.simout2struct(t,x,p);
                
                output_gray{i,j,k,l} = out;
            end
        end
    end
end

% Color.

output_color = {};

m = model_fedbatch();

for i = 1:size(input_color,1)
    for j = 1:size(input_color,2)
        for k = 1:size(input_color,3)
            for l = 1:size(input_color,4)
                % Default simulation parameters.
                p = m.p;
                
                % Change the parameters of the gene expression space.
                p.cell__p_A__k_u   = input_color{i,j,k,l}.p.cell__p_A__k_u;
                p.cell__p_A__k_b   = input_color{i,j,k,l}.p.cell__p_A__k_b;
                p.cell__p_A__omega = input_color{i,j,k,l}.p.cell__p_A__omega;
                
                % Default initial conditions.
                x0 = m.x0;
                
                % Use initial conditions of previous simulation.
                x0(1)  = input_color{i,j,k,l}.x0.p_r__m;
                x0(2)  = input_color{i,j,k,l}.x0.p_nr__m;
                x0(3)  = input_color{i,j,k,l}.x0.mu;
                x0(4)  = input_color{i,j,k,l}.x0.r;
                x0(5)  = input_color{i,j,k,l}.x0.p_A__m;
                x0(10) = input_color{i,j,k,l}.p.bio__s;
                
                % Solver options.
                opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
                opt = odeset(opt,'Mass',m.M);
                opt = odeset(opt,'Events',@(t,y) eventSubstrateDepletion(t,y,p,m));
                
                % Simulation time span.
                tspan = [0 1e9];
                
                % Simulate.
                [t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x0,opt);
                out = m.simout2struct(t,x,p);
                
                output_color{i,j,k,l} = out;
            end
        end
    end
end

%% Step 5: Calculate TRY.

% Gray.

for i = 1:size(input_gray,1)
    for j = 1:size(input_gray,2)
        for k = 1:size(input_gray,3)
            for l = 1:size(input_gray,4)
            V_0     = output_gray{i,j,k,l}.bio__V(1); % (L)
            V_f     = output_gray{i,j,k,l}.bio__V(end); % (L)
            Vout_f  = output_gray{i,j,k,l}.bio__V_out(end); % (L)
            Vfeed_f = output_gray{i,j,k,l}.bio__V_feed(end); % (L)
            
            n_f     = output_gray{i,j,k,l}.bio__N(end); % (teracells/L)
            
            mA_f    = output_gray{i,j,k,l}.cell__p_A__m(end); % (fg/cell)
            MA_f    = output_gray{i,j,k,l}.M_A(end); % (g)
            
            t_f     = output_gray{i,j,k,l}.t(end)/60; % (h)
            
            s_0     = output_gray{i,j,k,l}.bio__s(1); % (g/L)
            s_f     = output_gray{i,j,k,l}.bio__s(end); % (g/L)
            s_feed  = output_gray{i,j,k,l}.bio__s_f(end); % (g/L)
            Sout_f  = output_gray{i,j,k,l}.bio__S(end); % (g/L)
            
            
            titer   = (V_f*n_f*mA_f*1e-3+MA_f)/(V_f+Vout_f);
            productivity = titer/t_f;
            yield   = titer*(V_f+Vout_f)/(s_0*V_0 - s_f*V_f + s_feed*Vfeed_f - Sout_f);
            
            output_gray{i,j,k,l}.titer = titer;
            output_gray{i,j,k,l}.productivity = productivity;
            output_gray{i,j,k,l}.yield = yield;
            end
        end
    end
end

% Color.

for i = 1:size(input_color,1)
    for j = 1:size(input_color,2)
        for k = 1:size(input_color,3)
            for l = 1:size(input_color,4)
                V_0     = output_color{i,j,k,l}.bio__V(1); % (L)
                V_f     = output_color{i,j,k,l}.bio__V(end); % (L)
                Vout_f  = output_color{i,j,k,l}.bio__V_out(end); % (L)
                Vfeed_f = output_color{i,j,k,l}.bio__V_feed(end); % (L)
                
                n_f     = output_color{i,j,k,l}.bio__N(end); % (teracells/L)
                
                mA_f    = output_color{i,j,k,l}.cell__p_A__m(end); % (fg/cell)
                MA_f    = output_color{i,j,k,l}.M_A(end); % (g)
                
                t_f     = output_color{i,j,k,l}.t(end)/60; % (h)
                
                s_0     = output_color{i,j,k,l}.bio__s(1); % (g/L)
                s_f     = output_color{i,j,k,l}.bio__s(end); % (g/L)
                s_feed  = output_color{i,j,k,l}.bio__s_f(end); % (g/L)
                Sout_f  = output_color{i,j,k,l}.bio__S(end); % (g/L)
                
                
                titer   = (V_f*n_f*mA_f*1e-3+MA_f)/(V_f+Vout_f);
                productivity = titer/t_f;
                yield   = titer*(V_f+Vout_f)/(s_0*V_0 - s_f*V_f + s_feed*Vfeed_f - Sout_f);
                
                output_color{i,j,k,l}.titer = titer;
                output_color{i,j,k,l}.productivity = productivity;
                output_color{i,j,k,l}.yield = yield;
            end
        end
    end
end

%% Step 6: Calculate TRY relative variation indices.

% Gray.

for i = 1:size(input_gray,1)
    for j = 1:size(input_gray,2)
        for k = 1:size(input_gray,3)
            s = [];
            titer = [];
            productivity = [];
            yield = [];
            
            for l = 1:size(input_gray,4)
                s(l) = input_gray{i,j,k,l}.p.bio__s;
                titer(l) = output_gray{i,j,k,l}.titer;
                productivity(l) = output_gray{i,j,k,l}.productivity;
                yield(l) = output_gray{i,j,k,l}.yield;
            end
            
            titer_norm = titer./titer(end);
            productivity_norm = productivity./productivity(end);
            yield_norm = yield./yield(end);
            
            titer_variation = trapz(s,abs(titer_norm-1));
            productivity_variation = trapz(s,abs(productivity_norm-1));
            yield_variation = trapz(s,abs(yield_norm-1));
            
            output_gray{i,j,k,l}.titer_variation = titer_variation;
            output_gray{i,j,k,l}.productivity_variation = productivity_variation;
            output_gray{i,j,k,l}.yield_variation = yield_variation;
        end
    end
end

% Color.

for i = 1:size(input_color,1)
    for j = 1:size(input_color,2)
        for k = 1:size(input_color,3)
            s = [];
            titer = [];
            productivity = [];
            yield = [];
            
            for l = 1:size(input_color,4)
                s(l) = input_color{i,j,k,l}.p.bio__s;
                titer(l) = output_color{i,j,k,l}.titer;
                productivity(l) = output_color{i,j,k,l}.productivity;
                yield(l) = output_color{i,j,k,l}.yield;
            end
            
            titer_norm = titer./titer(end);
            productivity_norm = productivity./productivity(end);
            yield_norm = yield./yield(end);
            
            titer_variation = trapz(s,abs(titer_norm-1));
            productivity_variation = trapz(s,abs(productivity_norm-1));
            yield_variation = trapz(s,abs(yield_norm-1));
            
            output_color{i,j,k,l}.titer_variation = titer_variation;
            output_color{i,j,k,l}.productivity_variation = productivity_variation;
            output_color{i,j,k,l}.yield_variation = yield_variation;
        end
    end
end

%% Step 7: Prepare the level diagram data.

% Gray.

ld_gray.titer = [];
ld_gray.productivity = [];
ld_gray.yield = [];
ld_gray.titer_variation = [];
ld_gray.productivity_variation = [];
ld_gray.yield_variation = [];
ld_gray.omega = [];
ld_gray.K_C0 = [];
ld_gray.mu = [];

l = size(input_gray,4);

for i = 1:size(input_gray,1)
    for j = 1:size(input_gray,2)
        for k = 1:size(input_gray,3)
            ld_gray.titer(end+1) = output_gray{i,j,k,l}.titer;
            ld_gray.productivity(end+1) = output_gray{i,j,k,l}.productivity;
            ld_gray.yield(end+1) = output_gray{i,j,k,l}.yield;
            
            ld_gray.titer_variation(end+1) = output_gray{i,j,k,l}.titer_variation;
            ld_gray.productivity_variation(end+1) = output_gray{i,j,k,l}.productivity_variation;
            ld_gray.yield_variation(end+1) = output_gray{i,j,k,l}.yield_variation;
            
            
            ld_gray.omega(end+1) = output_gray{i,j,k,l}.cell__p_A__omega(end)*output_gray{i,j,k,l}.cell__p_A__N(end);
            ld_gray.K_C0(end+1) = mean(output_gray{i,j,k,l}.cell__p_A__K_C0(:));
            ld_gray.mu(end+1) = mean(output_gray{i,j,k,l}.cell__mu(:));
        end
    end
end

ld_gray.titer_norm = ld_gray.titer./max(ld_gray.titer);
ld_gray.productivity_norm = ld_gray.productivity./max(ld_gray.productivity);
ld_gray.yield_norm = ld_gray.yield./max(ld_gray.yield);

% Color.

ld_color.titer = [];
ld_color.productivity = [];
ld_color.yield = [];
ld_color.titer_variation = [];
ld_color.productivity_variation = [];
ld_color.yield_variation = [];
ld_color.omega = [];
ld_color.K_C0 = [];
ld_color.mu = [];
ld_color.color = {};
ld_color.size = [];

l = size(input_color,4);

for i = 1:size(input_color,1)
    for j = 1:size(input_color,2)
        for k = 1:size(input_color,3)
            ld_color.titer(end+1) = output_color{i,j,k,l}.titer;
            ld_color.productivity(end+1) = output_color{i,j,k,l}.productivity;
            ld_color.yield(end+1) = output_color{i,j,k,l}.yield;
            
            ld_color.titer_variation(end+1) = output_color{i,j,k,l}.titer_variation;
            ld_color.productivity_variation(end+1) = output_color{i,j,k,l}.productivity_variation;
            ld_color.yield_variation(end+1) = output_color{i,j,k,l}.yield_variation;
            
            ld_color.omega(end+1) = output_color{i,j,k,l}.cell__p_A__omega(end)*output_color{i,j,k,l}.cell__p_A__N(end);
            ld_color.K_C0(end+1) = mean(output_color{i,j,k,l}.cell__p_A__K_C0(1));
            ld_color.mu(end+1) = mean(output_color{i,j,k,l}.cell__mu(:));
            
            % Assign color based on omega value.
            if mod(k-1, 3) == 0
                ld_color.color{end+1} = [0, 0.4470, 0.7410]; % blue
            elseif mod(k-1, 3) == 1
                ld_color.color{end+1} = [0.9290, 0.6940, 0.1250]; %yellow
            else
                ld_color.color{end+1} = [0.6350, 0.0780, 0.1840]; %red
            end
            
            % Assign size based on kb value.
            if mod(j-1, 3) == 0
                ld_color.size(end+1) = 9;
            elseif mod(j-1, 3) == 1
                ld_color.size(end+1) = 15; 
            else
                ld_color.size(end+1) = 20;
            end
            
            
        end
    end
end

ld_color.titer_norm = ld_color.titer./max(ld_color.titer);
ld_color.productivity_norm = ld_color.productivity./max(ld_color.productivity);
ld_color.yield_norm = ld_color.yield./max(ld_color.yield);

%% Step 8: Plot results.

close all;

fig = figure('units','centimeters','position',[0 0 14.5 13]);

subplot(3,3,1);

hold on;
plot(ld_gray.titer, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.titer(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
ylabel('Growth rate [min$^{-1}$]','interpreter','latex','FontSize',9);
xlabel('Titer [g $\cdot$ L$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 80]);

subplot(3,3,2);

hold on;
plot(ld_gray.productivity, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.productivity(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Productivity [g $\cdot$ L$^{-1}$ $\cdot$ h$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);

subplot(3,3,3);

hold on;
plot(ld_gray.yield, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.yield(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Yield [g $\cdot$ g$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);


subplot(3,3,4);

hold on;
plot(ld_gray.titer_variation, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.titer_variation(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
ylabel('Growth rate [min$^{-1}$]','interpreter','latex','FontSize',9);
xlabel('Titer variation [adim]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 3]);

subplot(3,3,5);

hold on;
plot(ld_gray.productivity_variation, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.productivity_variation(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Productivity variation [adim]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 3]);

subplot(3,3,6);

hold on;
plot(ld_gray.yield_variation, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.yield_variation(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Yield variation [adim]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 3]);

subplot(3,3,7);

hold on;
plot(ld_gray.omega, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.omega(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
ylabel('Growth rate [min$^{-1}$]','interpreter','latex','FontSize',9);
xlabel('$N_A \omega_A$ [molec $\cdot$ min$^{-1}$ $\cdot$ cell$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);

subplot(3,3,8);

hold on;
plot(ld_gray.K_C0, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.K_C0(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('$K^A_{C_0}(s_n)$ [cell $\cdot$ molec$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 0.4]);

subplot(3,3,9);

hold on;
plot(ld_gray.mu, ld_gray.mu, '.', 'Color', [0.8 0.8 0.9]);
for i = length(ld_color.mu):-1:1
    plot(ld_color.mu(i), ld_color.mu(i), '.', 'Color', ld_color.color{i},'MarkerSize',ld_color.size(i));
end

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Mean growth rate [min$^{-1}$]','interpreter','latex','FontSize',9);
xlim([0 0.03]);
ylim([0 0.025]);

annotation('textbox', [0.05, 0.98, 0, 0], 'string', 'A', 'FontWeight', 'Bold','FontSize',11);
annotation('textbox', [0.05, 0.68, 0, 0], 'string', 'B', 'FontWeight', 'Bold','FontSize',11);
annotation('textbox', [0.05, 0.38, 0, 0], 'string', 'C', 'FontWeight', 'Bold','FontSize',11);

print(fig,'./figs/figure_3.eps','-depsc');
