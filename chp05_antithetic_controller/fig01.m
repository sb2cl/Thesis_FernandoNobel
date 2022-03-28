function [] = fig01(p_base,t_inc)
%% Parameters specific to this simulation.

p_wild = p_base;
p_wild.cell__A__omega_max = 0;
p_wild.cell__B__omega_max = 0;

p_open = p_base;
p_open.cell____openloop = 1;

p_close = p_base;
p_close.cell____openloop = 0;

%% Simulate.

out_wild  = simulate_host(p_wild, t_inc);
out_open  = simulate_host(p_open, t_inc);
out_close = simulate_host(p_close, t_inc);


%% Plot result.

lineWidth = 1.2;
colors;

fig = figure('units','centimeters','position',[0,13,13.5,9]);

subplot(2,2,1);

hold on; 
% grid on;

shadePerturbation(t_inc);
plot(out_open.t, out_open.cell__A__m,'LineWidth',lineWidth);
plot(out_close.t, out_close.cell__A__m,'LineWidth',lineWidth);
plot(out_close.t, out_close.cell__ref,'k:','LineWidth',lineWidth);
% plot(out.t, ones(size(out.t))*p_base.cell__x1__omega_max/p_base.cell__x2__omega_max,'k:','LineWidth',lineWidth);
ylim([0 12]);

set(gca,'Layer','top','FontSize',6);
ylabel('[fg $\cdot$ cell$^{-1}$]','interpreter','latex','FontSize',9);
title('\textbf{Output (Protein A)}','interpreter','latex','FontSize',9);
xlabel('Time [h]','interpreter','latex','FontSize',9);
leg = legend('open-loop', 'closed-loop', 'Location', 'southeast','interpreter','latex');
leg.ItemTokenSize = [10,10];

subplot(2,2,2);
hold on; 
% grid on;

shadePerturbation(t_inc);
plot(out_open.t, out_open.cell__B__m,'LineWidth',lineWidth);
plot(out_close.t, out_close.cell__B__m,'LineWidth',lineWidth);
ylim([0 50]);

set(gca,'Layer','top','FontSize',6);
ylabel('[fg $\cdot$ cell$^{-1}$]','interpreter','latex','FontSize',9);
title('\textbf{Perturbation (Protein B)}','interpreter','latex','FontSize',9);
xlabel('Time [h]','interpreter','latex','FontSize',9);

subplot(2,2,3);

hold on; 
% grid on;

shadePerturbation(t_inc);
plot(out_close.t, out_close.cell__x1__m,'LineWidth',lineWidth,'color',myColors.orange);
ylim([0 4]);

set(gca,'Layer','top','FontSize',6);
ylabel('[fg $\cdot$ cell$^{-1}$]','interpreter','latex','FontSize',9);
title('\textbf{Control action (Sigma factor)}','interpreter','latex','FontSize',9);
xlabel('Time [h]','interpreter','latex','FontSize',9);


subplot(2,2,4);

hold on; 
% grid on;

shadePerturbation(t_inc);
plot(out_open.t, out_open.cell__mu,'LineWidth',lineWidth);
plot(out_close.t, out_close.cell__mu,'LineWidth',lineWidth);
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
print(fig,'./figs/fig01.eps','-depsc');
end