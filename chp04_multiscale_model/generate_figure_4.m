close all;
clear all;

%% Step 1: Set of parameters for low, medium and high burden.

% all
kb = 15;
% s = [3.6];
% s = [0.1 1 2 3.6];
s  = [0.1 0.2 0.4 0.6 0.8 1 1.5 2 2.5 3 3.6];

% high burden
ku_high = [6 20 135];
productivity_desired_high = 1;

% low burden
ku_low = [6 20 135];
productivity_desired_low = 0.1;

input = {};

% High-burden.
for i = 1:3
    for j = 1:length(s)
        input{i,j}.p.cell__p_A__k_u = ku_high(i);
        input{i,j}.p.cell__p_A__k_b = kb;
        input{i,j}.p.bio__s = s(j);
        input{i,j}.productivity_desired = productivity_desired_high;
    end
end

% Medium-burden.
for i = 4:6
    for j = 1:length(s)
        input{i,j}.p.cell__p_A__k_u = ku_high(i-3);
        input{i,j}.p.cell__p_A__k_b = kb;
        input{i,j}.p.bio__s = s(j);
        input{i,j}.productivity_desired = productivity_desired_high;
    end
end

% Low-burden.
for i = 7:9
    for j = 1:length(s)
        input{i,j}.p.cell__p_A__k_u = ku_low(i-6);
        input{i,j}.p.cell__p_A__k_b = kb;
        input{i,j}.p.bio__s = s(j);
        input{i,j}.productivity_desired = productivity_desired_low;
    end
end

%% Step 2: Find the value of omega that ahieves the desired productivity.

omega = [];
omega(1) = 124.8;
omega(2) = 156.55;
omega(3) = 417.40;
omega(4) = 21.68;
omega(5) = 26.27;
omega(6) = 63.97;
omega(7) = 1.3763;
omega(8) = 1.507; % 2.0665;
omega(9) = 2.5795;

for i = 1:size(input,1)
    for j = 1:size(input,2)
        input{i,j}.p.cell__p_A__omega = omega(i);
    end
end

%% Step 3: Calculate initial conditions for the set of parameters.

% Init model.
m = model_initialize();

% Solver options.
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
opt = odeset(opt,'Mass',m.M);

% Simulation time span.
tspan = [m.opts.t_init m.opts.t_end];

for i = 1:size(input,1)
    for j = 1:size(input,2)
        % Default simulation parameters.
        p = m.p;
        
        % Change the parameters of the gene expression space.
        p.cell__p_A__k_u   = input{i,j}.p.cell__p_A__k_u;
        p.cell__p_A__k_b   = input{i,j}.p.cell__p_A__k_b;
        p.cell__p_A__omega = input{i,j}.p.cell__p_A__omega;
        p.s                = input{i,j}.p.bio__s;
        
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
        input{i,j}.x0 = x0;
        
    end
end

%% Step 4: Simulate the set with model_fedbatch.

output = {};

m = model_fedbatch();

for i = 1:size(input,1)
    for j = 1:size(input,2)
        % Default simulation parameters.
        p = m.p;
        
        % Change the parameters of the gene expression space.
        p.cell__p_A__k_u   = input{i,j}.p.cell__p_A__k_u;
        p.cell__p_A__k_b   = input{i,j}.p.cell__p_A__k_b;
        p.cell__p_A__omega = input{i,j}.p.cell__p_A__omega;
        
        % Default initial conditions.
        x0 = m.x0;
        
        % Use initial conditions of previous simulation.
        x0(1)  = input{i,j}.x0.p_r__m;
        x0(2)  = input{i,j}.x0.p_nr__m;
        x0(3)  = input{i,j}.x0.mu;
        x0(4)  = input{i,j}.x0.r;
        x0(5)  = input{i,j}.x0.p_A__m;
        x0(10) = input{i,j}.p.bio__s;
        
        % Solver options.
        opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
        opt = odeset(opt,'Mass',m.M);
        opt = odeset(opt,'Events',@(t,y) eventSubstrateDepletion(t,y,p,m));
        
        % Simulation time span.
        tspan = [0 1e9];
        
        % Simulate.
        [t,x] = ode15s(@(t,x) m.ode(t,x,p),tspan,x0,opt);
        out = m.simout2struct(t,x,p);
        
        output{i,j} = out;
    end
end

%% Step 5: Calculate TRY.

for i = 1:size(input,1)
    for j = 1:size(input,2)
            V_0     = output{i,j}.bio__V(1); % (L)
            V_f     = output{i,j}.bio__V(end); % (L)
            Vout_f  = output{i,j}.bio__V_out(end); % (L)
            Vfeed_f = output{i,j}.bio__V_feed(end); % (L)
            
            n_f     = output{i,j}.bio__N(end); % (teracells/L)
            
            mA_f    = output{i,j}.cell__p_A__m(end); % (fg/cell)
            MA_f    = output{i,j}.M_A(end); % (g)
            
            t_f     = output{i,j}.t(end)/60; % (h)
            
            s_0     = output{i,j}.bio__s(1); % (g/L)
            s_f     = output{i,j}.bio__s(end); % (g/L)
            s_feed  = output{i,j}.bio__s_f(end); % (g/L)
            Sout_f  = output{i,j}.bio__S(end); % (g/L)
            
            
            titer   = (V_f*n_f*mA_f*1e-3+MA_f)/(V_f+Vout_f);
            productivity = titer/t_f;
            yield   = titer*(V_f+Vout_f)/(s_0*V_0 - s_f*V_f + s_feed*Vfeed_f - Sout_f);
            
            output{i,j}.titer = titer;
            output{i,j}.productivity = productivity;
            output{i,j}.yield = yield;
    end
end

%% Step 6: Prepare the level diagram data.

ld = {};

for i = 1:size(input,1)
    ld{i}.titer = [];
    ld{i}.productivity = [];
    ld{i}.yield = [];
    ld{i}.mu = [];
    ld{i}.substrate = [];
    
    for j = 1:size(input,2)
            ld{i}.titer(end+1) = output{i,j}.titer;
            ld{i}.productivity(end+1) = output{i,j}.productivity;
            ld{i}.yield(end+1) = output{i,j}.yield;
            ld{i}.mu(end+1) = mean(output{i,j}.cell__mu);
            ld{i}.substrate(end+1) = output{i,j}.bio__s(1);
    end
end

%% Step 7: Plot results.

lineWidth = 1.2;

close all;

fig = figure('units','centimeters','position',[0 0 14.5 13]);

subplot(4,2,1);
hold on;
plot(ld{1}.substrate,ld{1}.titer,'Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{2}.substrate,ld{2}.titer,'--','Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{3}.substrate,ld{3}.titer,':','Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{4}.substrate,ld{4}.titer,'Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);
plot(ld{5}.substrate,ld{5}.titer,'--','Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);
plot(ld{6}.substrate,ld{6}.titer,':','Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);

grid on;
ax = gca;
ax.FontSize = 6;
ylabel('Titer [g $\cdot$ L$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylim([0 60]);

subplot(4,2,3);
hold on;
plot(ld{1}.substrate,ld{1}.productivity,'Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{2}.substrate,ld{2}.productivity,'--','Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{3}.substrate,ld{3}.productivity,':','Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{4}.substrate,ld{4}.productivity,'Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);
plot(ld{5}.substrate,ld{5}.productivity,'--','Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);
plot(ld{6}.substrate,ld{6}.productivity,':','Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);

grid on;
ax = gca;
ax.FontSize = 6;
ylabel('Productivity [g $\cdot$ L$^{-1}$ $\cdot$ h$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylim([0 1.5]);

subplot(4,2,5);
hold on;
plot(ld{1}.substrate,ld{1}.yield,'Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{2}.substrate,ld{2}.yield,'--','Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{3}.substrate,ld{3}.yield,':','Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{4}.substrate,ld{4}.yield,'Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);
plot(ld{5}.substrate,ld{5}.yield,'--','Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);
plot(ld{6}.substrate,ld{6}.yield,':','Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);

grid on;
ax = gca;
ax.FontSize = 6;
ylabel('Yield [g $\cdot$ g$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylim([0 +0.4]);

subplot(4,2,7);
hold on;
plot(ld{1}.substrate,ld{1}.mu,'Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{2}.substrate,ld{2}.mu,'--','Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{3}.substrate,ld{3}.mu,':','Color', [0.8500 0.3250 0.0980],'LineWidth',lineWidth);
plot(ld{4}.substrate,ld{4}.mu,'Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);
plot(ld{5}.substrate,ld{5}.mu,'--','Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);
plot(ld{6}.substrate,ld{6}.mu,':','Color', [0.9290 0.6940 0.1250],'LineWidth',lineWidth);

grid on;
ax = gca;
ax.FontSize = 6;
xlabel('Substrate [g $\cdot$ L$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylabel('Growth rate [min$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylim([0 0.015]);

subplot(4,2,2);
hold on;
plot(ld{7}.substrate,ld{7}.titer,'Color', [0 0.4470 0.7410],'LineWidth',lineWidth);
plot(ld{8}.substrate,ld{8}.titer,'--','Color', [0 0.4470 0.7410],'LineWidth',lineWidth);
plot(ld{9}.substrate,ld{9}.titer,':','Color', [0 0.4470 0.7410],'LineWidth',lineWidth);

grid on;
ax = gca;
ax.FontSize = 6;
ylabel('Titer [g $\cdot$ L$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylim([0 3]);

subplot(4,2,4);
hold on;
plot(ld{7}.substrate,ld{7}.productivity,'Color', [0 0.4470 0.7410],'LineWidth',lineWidth);
plot(ld{8}.substrate,ld{8}.productivity,'--','Color', [0 0.4470 0.7410],'LineWidth',lineWidth);
plot(ld{9}.substrate,ld{9}.productivity,':','Color', [0 0.4470 0.7410],'LineWidth',lineWidth);

grid on;
ax = gca;
ax.FontSize = 6;
ylabel('Productivity [g $\cdot$ L$^{-1}$ $\cdot$ h$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylim([0 0.15]);

subplot(4,2,6);
hold on;
plot(ld{7}.substrate,ld{7}.yield,'Color', [0 0.4470 0.7410],'LineWidth',lineWidth);
plot(ld{8}.substrate,ld{8}.yield,'--','Color', [0 0.4470 0.7410],'LineWidth',lineWidth);
plot(ld{9}.substrate,ld{9}.yield,':','Color', [0 0.4470 0.7410],'LineWidth',lineWidth);

grid on;
ax = gca;
ax.FontSize = 6;
ylabel('Yield [g $\cdot$ g$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylim([0 0.02]);

subplot(4,2,8);
hold on;
plot(ld{7}.substrate,ld{7}.mu,'Color', [0 0.4470 0.7410],'LineWidth',lineWidth);
plot(ld{8}.substrate,ld{8}.mu,'--','Color', [0 0.4470 0.7410],'LineWidth',lineWidth);
plot(ld{9}.substrate,ld{9}.mu,':','Color', [0 0.4470 0.7410],'LineWidth',lineWidth);

grid on;
ax = gca;
ax.FontSize = 6;
ylabel('Growth rate [min$^{-1}$]', 'interpreter', 'latex','FontSize',8);
xlabel('Substrate [g $\cdot$ L$^{-1}$]', 'interpreter', 'latex','FontSize',8);
ylim([0 0.02]);

annotation('textbox', [0.05, 0.98, 0, 0], 'string', 'A', 'FontWeight', 'Bold','FontSize',11);
annotation('textbox', [0.50, 0.98, 0, 0], 'string', 'B', 'FontWeight', 'Bold','FontSize',11);

print(fig,'./figs/figure_4.eps','-depsc');
