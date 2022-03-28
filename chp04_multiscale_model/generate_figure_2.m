clear all;
close all;

%% Step 1: Generate a 3-D array of values ku, kb, and N*omega.

% Limit values of gene expression space.
lim_ku    = [6 135];
lim_kb    = [3 15];
lim_omega = [100 350];

% Number of values picked between the limits.
n_ku    = 4;
n_kb    = 4;
n_omega = 4;

% Generate array with the values choosen to form the gene expression space.
ku    = arrayMinMaxN(lim_ku, n_ku);
kb    = arrayMinMaxN(lim_kb, n_kb);
omega = arrayMinMaxN(lim_omega, n_omega);

% Gene expression space.
input = {};

% Fill the gene expression space.
for i_ku = 1:n_ku
    for i_kb = 1:n_kb
        for i_omega = 1:n_omega
            input{i_ku,i_kb,i_omega}.p.cell__p_A__k_u = ku(i_ku);
            input{i_ku,i_kb,i_omega}.p.cell__p_A__k_b = kb(i_kb);
            input{i_ku,i_kb,i_omega}.p.cell__p_A__omega = omega(i_omega);
        end
    end
end

%% Step 2: Calculate initial conditions with the unlimited substrate model.

% Init model.
m = model_initialize();

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [m.opts.t_init m.opts.t_end];

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            % Default simulation parameters.
            p = m.p;
            
            % Change the parameters of the gene expression space.
            p.cell__p_A__k_u   = input{i,j,k}.p.cell__p_A__k_u;
            p.cell__p_A__k_b   = input{i,j,k}.p.cell__p_A__k_b;
            p.cell__p_A__omega = input{i,j,k}.p.cell__p_A__omega;
    
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
            input{i,j,k}.x0 = x0;
        end
    end
end

%% Step 3: Simulate batch.

output = {};

m = model_batch();

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            % Default simulation parameters.
            p = m.p;
            
            % Change the parameters of the gene expression space.
            p.cell__p_A__k_u   = input{i,j,k}.p.cell__p_A__k_u;
            p.cell__p_A__k_b   = input{i,j,k}.p.cell__p_A__k_b;
            p.cell__p_A__omega = input{i,j,k}.p.cell__p_A__omega;
            
            % Default initial conditions.
            x0 = m.x0;
            
            % Use initial conditions of previous simulation.
            x0(1) = input{i,j,k}.x0.p_r__m;
            x0(2) = input{i,j,k}.x0.p_nr__m;
            x0(3) = input{i,j,k}.x0.mu;
            x0(4) = input{i,j,k}.x0.r;
            x0(5) = input{i,j,k}.x0.p_A__m;
            
            % Solver options.
            opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
            opt = odeset(opt,'Mass',m.M);
            opt = odeset(opt,'Events',@(t,y) eventSubstrateDepletion(t,y,p,m));
            
            % Simulation time span.
            tspan = [0 1e9];

            % Simulate.
            [t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x0,opt);
            out = m.simout2struct(t,x,p);  
            
            output{i,j,k}.batch = out;
        end
    end
end

% for i = 1:n_ku
%     for j = 1:n_kb
%         for k = 1:n_omega
% 
% out = output{i,j,k}.batch;
% 
% figure(1);
% 
% hold on;
% plot(out.t, out.bio__s);
% plot(out.t(end), out.bio__s(end),'o');
% 
%         end
%     end
% end

%% Step 4: Simulate fedbacth.

m = model_fedbatch();

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            % Default simulation parameters.
            p = m.p;
            
            % Change the parameters of the gene expression space.
            p.cell__p_A__k_u   = input{i,j,k}.p.cell__p_A__k_u;
            p.cell__p_A__k_b   = input{i,j,k}.p.cell__p_A__k_b;
            p.cell__p_A__omega = input{i,j,k}.p.cell__p_A__omega;
            
            % Default initial conditions.
            x0 = m.x0;
            
            % Use initial conditions of previous simulation.
            x0(1) = input{i,j,k}.x0.p_r__m;
            x0(2) = input{i,j,k}.x0.p_nr__m;
            x0(3) = input{i,j,k}.x0.mu;
            x0(4) = input{i,j,k}.x0.r;
            x0(5) = input{i,j,k}.x0.p_A__m;
            
            % Solver options.
            opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
            opt = odeset(opt,'Mass',m.M);
            opt = odeset(opt,'Events',@(t,y) eventSubstrateDepletion(t,y,p,m));
            
            % Simulation time span.
            tspan = [0 1e9];

            % Simulate.
            [t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x0,opt);
            out = m.simout2struct(t,x,p);  
            
            output{i,j,k}.fedbatch = out;
        end
    end
end

% for i = 1:n_ku
%     for j = 1:n_kb
%         for k = 1:n_omega
% 
% out = output{i,j,k}.fedbatch;
% 
% figure(1);
% 
% hold on;
% plot(out.t, out.bio__s);
% plot(out.t(end), out.bio__s(end),'o');
% 
%         end
%     end
% end

%% Step 5: Simulate continous.

m = model_continous();

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            % Default simulation parameters.
            p = m.p;
            
            % Change the parameters of the gene expression space.
            p.cell__p_A__k_u   = input{i,j,k}.p.cell__p_A__k_u;
            p.cell__p_A__k_b   = input{i,j,k}.p.cell__p_A__k_b;
            p.cell__p_A__omega = input{i,j,k}.p.cell__p_A__omega;
            
            % Default initial conditions.
            x0 = m.x0;
            
            % Use initial conditions of previous simulation.
            x0(1) = input{i,j,k}.x0.p_r__m;
            x0(2) = input{i,j,k}.x0.p_nr__m;
            x0(3) = input{i,j,k}.x0.mu;
            x0(4) = input{i,j,k}.x0.r;
            x0(5) = input{i,j,k}.x0.p_A__m;
            
            % Solver options.
            opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
            opt = odeset(opt,'Mass',m.M);
            opt = odeset(opt,'Events',@(t,y) eventFeedDepletion(t,y,p,m));
            
            % Simulation time span.
            tspan = [0 1e9];

            % Simulate.
            [t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x0,opt);
            out = m.simout2struct(t,x,p);  
            
            output{i,j,k}.continous = out;
        end
    end
end

% for i = 1:n_ku
%     for j = 1:n_kb
%         for k = 1:n_omega
% 
% out = output{i,j,k}.continous;
% 
% figure(1);
% 
% hold on;
% plot(out.t, out.bio__s);
% plot(out.t(end), out.bio__s(end),'o');
% 
% figure(2);
% 
% hold on;
% plot(out.t, out.bio__x);
% plot(out.t(end), out.bio__x(end),'o');
% 
%         end
%     end
% end

%% Step 6: Calculate TRY.

% Batch

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            V_0     = output{i,j,k}.batch.bio__V(1); % (L)
            V_f     = output{i,j,k}.batch.bio__V(end); % (L)
            Vout_f  = output{i,j,k}.batch.bio__V_out(end); % (L)
            Vfeed_f = output{i,j,k}.batch.bio__V_feed(end); % (L)
            
            n_f     = output{i,j,k}.batch.bio__N(end); % (teracells/L)
            
            mA_f    = output{i,j,k}.batch.cell__p_A__m(end); % (fg/cell)
            MA_f    = output{i,j,k}.batch.M_A(end); % (g)
            
            t_f     = output{i,j,k}.batch.t(end)/60; % (h)
            
            s_0     = output{i,j,k}.batch.bio__s(1); % (g/L)
            s_f     = output{i,j,k}.batch.bio__s(end); % (g/L)
            s_feed  = output{i,j,k}.batch.bio__s_f(end); % (g/L)
            Sout_f  = output{i,j,k}.batch.bio__S(end); % (g/L)
            
            
            titer   = (V_f*n_f*mA_f*1e-3+MA_f)/(V_f+Vout_f);
            productivity = titer/t_f;
            yield   = titer*(V_f+Vout_f)/(s_0*V_0 - s_f*V_f + s_feed*Vfeed_f - Sout_f);
            
            output{i,j,k}.batch.titer = titer;
            output{i,j,k}.batch.productivity = productivity;
            output{i,j,k}.batch.yield = yield;
        end
    end
end

% Fedbatch

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            V_0     = output{i,j,k}.fedbatch.bio__V(1); % (L)
            V_f     = output{i,j,k}.fedbatch.bio__V(end); % (L)
            Vout_f  = output{i,j,k}.fedbatch.bio__V_out(end); % (L)
            Vfeed_f = output{i,j,k}.fedbatch.bio__V_feed(end); % (L)
            
            n_f     = output{i,j,k}.fedbatch.bio__N(end); % (teracells/L)
            
            mA_f    = output{i,j,k}.fedbatch.cell__p_A__m(end); % (fg/cell)
            MA_f    = output{i,j,k}.fedbatch.M_A(end); % (g)
            
            t_f     = output{i,j,k}.fedbatch.t(end)/60; % (h)
            
            s_0     = output{i,j,k}.fedbatch.bio__s(1); % (g/L)
            s_f     = output{i,j,k}.fedbatch.bio__s(end); % (g/L)
            s_feed  = output{i,j,k}.fedbatch.bio__s_f(end); % (g/L)
            Sout_f  = output{i,j,k}.fedbatch.bio__S(end); % (g/L)
            
            
            titer   = (V_f*n_f*mA_f*1e-3+MA_f)/(V_f+Vout_f);
            productivity = titer/t_f;
            yield   = titer*(V_f+Vout_f)/(s_0*V_0 - s_f*V_f + s_feed*Vfeed_f - Sout_f);
            
            output{i,j,k}.fedbatch.titer = titer;
            output{i,j,k}.fedbatch.productivity = productivity;
            output{i,j,k}.fedbatch.yield = yield;
        end
    end
end

% Continous

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            V_0     = output{i,j,k}.continous.bio__V(1); % (L)
            V_f     = output{i,j,k}.continous.bio__V(end); % (L)
            Vout_f  = output{i,j,k}.continous.bio__V_out(end); % (L)
            Vfeed_f = output{i,j,k}.continous.bio__V_feed(end); % (L)
            
            n_f     = output{i,j,k}.continous.bio__N(end); % (teracells/L)
            
            mA_f    = output{i,j,k}.continous.cell__p_A__m(end); % (fg/cell)
            MA_f    = output{i,j,k}.continous.M_A(end); % (g)
            
            t_f     = output{i,j,k}.continous.t(end)/60; % (h)
            
            s_0     = output{i,j,k}.continous.bio__s(1); % (g/L)
            s_f     = output{i,j,k}.continous.bio__s(end); % (g/L)
            s_feed  = output{i,j,k}.continous.bio__s_f(end); % (g/L)
            Sout_f  = output{i,j,k}.continous.bio__S(end); % (g/L)
            
            
            titer   = (V_f*n_f*mA_f*1e-3+MA_f)/(V_f+Vout_f);
            productivity = titer/t_f;
            yield   = titer*(V_f+Vout_f)/(s_0*V_0 - s_f*V_f + s_feed*Vfeed_f - Sout_f);
            
            output{i,j,k}.continous.titer = titer;
            output{i,j,k}.continous.productivity = productivity;
            output{i,j,k}.continous.yield = yield;
        end
    end
end

%% Step 7: Prepare the level diagram data.

% Batch

ld.batch.titer = [];
ld.batch.productivity = [];
ld.batch.yield = [];
ld.batch.omega = [];
ld.batch.K_C0 = [];
ld.batch.mu = [];

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            ld.batch.titer(end+1) = output{i,j,k}.batch.titer;
            ld.batch.productivity(end+1) = output{i,j,k}.batch.productivity;
            ld.batch.yield(end+1) = output{i,j,k}.batch.yield;
            ld.batch.omega(end+1) = output{i,j,k}.batch.cell__p_A__omega(end)*output{i,j,k}.batch.cell__p_A__N(end);
            ld.batch.K_C0(end+1) = mean(output{i,j,k}.batch.cell__p_A__K_C0(:));
            ld.batch.mu(end+1) = mean(output{i,j,k}.batch.cell__mu(:));
        end
    end
end

ld.batch.titer_norm = ld.batch.titer./max(ld.batch.titer);
ld.batch.productivity_norm = ld.batch.productivity./max(ld.batch.productivity);
ld.batch.yield_norm = ld.batch.yield./max(ld.batch.yield);

% Fedbatch

ld.fedbatch.titer = [];
ld.fedbatch.productivity = [];
ld.fedbatch.yield = [];
ld.fedbatch.omega = [];
ld.fedbatch.K_C0 = [];
ld.fedbatch.mu = [];

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            ld.fedbatch.titer(end+1) = output{i,j,k}.fedbatch.titer;
            ld.fedbatch.productivity(end+1) = output{i,j,k}.fedbatch.productivity;
            ld.fedbatch.yield(end+1) = output{i,j,k}.fedbatch.yield;
            ld.fedbatch.omega(end+1) = output{i,j,k}.fedbatch.cell__p_A__omega(end)*output{i,j,k}.fedbatch.cell__p_A__N(end);
            ld.fedbatch.K_C0(end+1) = mean(output{i,j,k}.fedbatch.cell__p_A__K_C0(:));
            ld.fedbatch.mu(end+1) = mean(output{i,j,k}.fedbatch.cell__mu(:));
        end
    end
end

ld.fedbatch.titer_norm = ld.fedbatch.titer./max(ld.fedbatch.titer);
ld.fedbatch.productivity_norm = ld.fedbatch.productivity./max(ld.fedbatch.productivity);
ld.fedbatch.yield_norm = ld.fedbatch.yield./max(ld.fedbatch.yield);

% Continous

ld.continous.titer = [];
ld.continous.productivity = [];
ld.continous.yield = [];
ld.continous.omega = [];
ld.continous.K_C0 = [];
ld.continous.mu = [];

for i = 1:n_ku
    for j = 1:n_kb
        for k = 1:n_omega
            ld.continous.titer(end+1) = output{i,j,k}.continous.titer;
            ld.continous.productivity(end+1) = output{i,j,k}.continous.productivity;
            ld.continous.yield(end+1) = output{i,j,k}.continous.yield;
            ld.continous.omega(end+1) = output{i,j,k}.continous.cell__p_A__omega(end)*output{i,j,k}.continous.cell__p_A__N(end);
            ld.continous.K_C0(end+1) = mean(output{i,j,k}.continous.cell__p_A__K_C0(1));
            ld.continous.mu(end+1) = mean(output{i,j,k}.continous.cell__mu(:));
        end
    end
end

ld.continous.titer_norm = ld.continous.titer./max(ld.continous.titer);
ld.continous.productivity_norm = ld.continous.productivity./max(ld.continous.productivity);
ld.continous.yield_norm = ld.continous.yield./max(ld.continous.yield);

%% Step 8: Plot results.

close all;

fig = figure('units','centimeters','position',[0 0 14.5 13]);

subplot(3,3,1);

hold on;
plot(ld.batch.titer, ld.batch.mu, '.');
plot(ld.fedbatch.titer, ld.fedbatch.mu, '.');
plot(ld.continous.titer, ld.continous.mu, '.');

grid on;
ax = gca;
ax.FontSize = 6; 
ylabel('Growth rate [min$^{-1}$]','interpreter','latex','FontSize',9);
xlabel('Titer [g $\cdot$ L$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 80]);


subplot(3,3,2);

hold on;
plot(ld.batch.productivity, ld.batch.mu, '.');
plot(ld.fedbatch.productivity, ld.fedbatch.mu, '.');
plot(ld.continous.productivity, ld.continous.mu, '.');

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Productivity [g $\cdot$ L$^{-1}$ $\cdot$ h$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 1.5]);

subplot(3,3,3);

hold on;
plot(ld.batch.yield, ld.batch.mu, '.');
plot(ld.fedbatch.yield, ld.fedbatch.mu, '.');
plot(ld.continous.yield, ld.continous.mu, '.');

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Yield [g $\cdot$ g$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 0.4]);

subplot(3,3,4);

hold on;
plot(ld.batch.titer_norm, ld.batch.mu, '.');
plot(ld.fedbatch.titer_norm, ld.fedbatch.mu, '.');
plot(ld.continous.titer_norm, ld.continous.mu, '.');

grid on;
ax = gca;
ax.FontSize = 6; 
ylabel('Growth rate [min$^{-1}$]','interpreter','latex','FontSize',9);
xlabel('Normalized titer [adim]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 +inf]);

subplot(3,3,5);

hold on;
plot(ld.batch.productivity_norm, ld.batch.mu, '.');
plot(ld.fedbatch.productivity_norm, ld.fedbatch.mu, '.');
plot(ld.continous.productivity_norm, ld.continous.mu, '.');

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Normalized productivity [adim]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 +inf]);

subplot(3,3,6);

hold on;
plot(ld.batch.yield_norm, ld.batch.mu, '.');
plot(ld.fedbatch.yield_norm, ld.fedbatch.mu, '.');
plot(ld.continous.yield_norm, ld.continous.mu, '.');

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Normalized yield [adim]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 +inf]);

subplot(3,3,7);

hold on;
plot(ld.batch.omega, ld.batch.mu, '.');
plot(ld.fedbatch.omega, ld.fedbatch.mu, '.');
plot(ld.continous.omega, ld.continous.mu, '.');

grid on;
ax = gca;
ax.FontSize = 6; 
ylabel('Growth rate [min$^{-1}$]','interpreter','latex','FontSize',9);
xlabel('$N_A \omega_A$ [molec $\cdot$ min$^{-1}$ $\cdot$ cell$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 400]);

subplot(3,3,8);

hold on;
plot(ld.batch.K_C0, ld.batch.mu, '.');
plot(ld.fedbatch.K_C0, ld.fedbatch.mu, '.');
plot(ld.continous.K_C0, ld.continous.mu, '.');

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('$K^A_{C_0}(s_n)e$ [cell $\cdot$ molec$^{-1}$]','interpreter','latex','FontSize',9);
ylim([0 0.025]);
xlim([0 0.4]);

subplot(3,3,9);

hold on;
plot(ld.batch.mu, ld.batch.mu, '.');
plot(ld.fedbatch.mu, ld.fedbatch.mu, '.');
plot(ld.continous.mu, ld.continous.mu, '.');
plot(0.0224167, 0.0224167, 'ks','MarkerFaceColor',[0 0 0]);

grid on;
ax = gca;
ax.FontSize = 6; 
xlabel('Mean growth rate [min$^{-1}$]','interpreter','latex','FontSize',9);
xlim([0 0.03]);
ylim([0 0.025]);

annotation('textbox', [0.05, 0.98, 0, 0], 'string', 'A', 'FontWeight', 'Bold','FontSize',11);
annotation('textbox', [0.05, 0.68, 0, 0], 'string', 'B', 'FontWeight', 'Bold','FontSize',11);
annotation('textbox', [0.05, 0.38, 0, 0], 'string', 'C', 'FontWeight', 'Bold','FontSize',11);

print(fig,'./figs/figure_2.eps','-depsc');
