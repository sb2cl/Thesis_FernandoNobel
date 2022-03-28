function [] = fig02(p_base,t_inc)
%% Parameters specific to this simulation.

K = [1 5 50];

for i = 1:length(K)
    p{i} = p_base;
    p{i}.cell__x1__omega_max = p_base.cell__x1__omega_max*K(i);
    p{i}.cell__x2__omega_max = p_base.cell__x2__omega_max*K(i);
    p{i}.cell____openloop = 0;
end

p_wild = p_base;
p_wild.cell__A__omega_max = 0;
p_wild.cell__B__omega_max = 0;

%% Simulate.

for i = 1:length(p)
    out{i} = simulate_host(p{i}, t_inc);
end

out_wild  = simulate_host(p_wild, t_inc);

%% Plot result.

lineWidth = 1.2;
colors;

colors_{1} = myColors.blue;
colors_{2} = myColors.yellow;
colors_{3} = myColors.orange;

fig = figure('units','centimeters','position',[13.5,13,13.5,9]);

subplot(2,2,1);

hold on; 
% grid on;

shadePerturbation(t_inc);
for i = 1:length(out)
    plot(out{i}.t, out{i}.cell__A__m,'LineWidth',lineWidth,'color',colors_{i});
end
plot(out{i}.t, out{i}.cell__ref,'k:','LineWidth',lineWidth);
ylim([0 6]);

set(gca,'Layer','top','FontSize',6);
ylabel('[fg $\cdot$ cell$^{-1}$]','interpreter','latex','FontSize',9);
title('\textbf{Output (Protein A)}','interpreter','latex','FontSize',9);
xlabel('Time [h]','interpreter','latex','FontSize',9);
leg = legend('base gain', 'high gain', 'very high gain', 'Location', 'southeast','interpreter','latex');
leg.ItemTokenSize = [10,10];

subplot(2,2,2);

hold on; 
% grid on;

shadePerturbation(t_inc);
for i = 1:length(out)
    plot(out{i}.t, out{i}.cell__B__m,'LineWidth',lineWidth,'color',colors_{i});
end
ylim([0 50]);

set(gca,'Layer','top','FontSize',6);
ylabel('[fg $\cdot$ cell$^{-1}$]','interpreter','latex','FontSize',9);
title('\textbf{Perturbation (Protein B)}','interpreter','latex','FontSize',9);
xlabel('Time [h]','interpreter','latex','FontSize',9);

subplot(2,2,3);

hold on; 
% grid on;

shadePerturbation(t_inc);
for i = 1:length(out)
    plot(out{i}.t, out{i}.cell__x1__m,'LineWidth',lineWidth,'color',colors_{i});
end
ylim([0 50]);

set(gca,'Layer','top','FontSize',6);
ylabel('[fg $\cdot$ cell$^{-1}$]','interpreter','latex','FontSize',9);
title('\textbf{Control action (Sigma factor)}','interpreter','latex','FontSize',9);
xlabel('Time [h]','interpreter','latex','FontSize',9);


subplot(2,2,4);

hold on; 
% grid on;

shadePerturbation(t_inc);
for i = 1:length(out)
    plot(out{i}.t, out{i}.cell__mu,'LineWidth',lineWidth,'color',colors_{i});
end
plot(out_wild.t, out_wild.cell__mu,'k:','LineWidth',lineWidth);
ylim([0 0.03]);

set(gca,'Layer','top','FontSize',6);
ylabel('[min$^{-1}$]','interpreter','latex','FontSize',9);
title('\textbf{Growth rate}','interpreter','latex','FontSize',9);
xlabel('Time [h]','interpreter','latex','FontSize',9);

annotation('textbox', [0.05, 0.98, 0, 0], 'string', 'A', 'FontWeight', 'Bold','FontSize',11);
annotation('textbox', [0.48, 0.98, 0, 0], 'string', 'B', 'FontWeight', 'Bold','FontSize',11);
annotation('textbox', [0.05, 0.48, 0, 0], 'string', 'C', 'FontWeight', 'Bold','FontSize',11);
annotation('textbox', [0.48, 0.48, 0, 0], 'string', 'D', 'FontWeight', 'Bold','FontSize',11);

set(gcf,'renderer','Painters');
print(fig,'./figs/fig02.eps','-depsc');
end