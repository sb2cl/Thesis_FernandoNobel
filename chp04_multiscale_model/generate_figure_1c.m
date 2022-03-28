s = [0.25 1 2 3.6];
x = [0.25 0.8 1.1 1];
S = [0.25:0.1:3.6];
X = interp1(s,x,S,'spline');
Y = ones(size(S));

close all;
fig = figure('units','centimeters','position',[0 0 5.5 4]);

hold on;

patch([S fliplr(S)], [X fliplr(Y)],[1 0.9 0.8],'LineStyle','none')
plot(S,X,'k','LineWidth',1.2);
plot([0 3.6],[1 1],'k:');
plot([3.6 3.6],[1 0],'k:');
plot([0.25 0.25],[0.99 0],'k:');
ylim([0 1.25]);

yticks([1]);
yticklabels({'$\mathcal{X}(s_n)$'});

xticks([0.25 3.6]);
xticklabels({'$s_{min}$','$s_n$'});
set(gca,'TickLabelInterpreter','latex');

xlabel('Substrate','interpreter','latex','FontSize',9);
print(fig,'./figs/figure_1c.eps','-depsc');


