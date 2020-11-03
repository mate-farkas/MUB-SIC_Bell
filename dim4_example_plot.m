load('dim4_example_data.mat')
plot(x, e, 'LineWidth', 1)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi / 2', '\pi', '3\pi / 2', '2\pi'})
xlim([0 2*pi])
xlabel('x','FontSize',14)
ylabel('\lambda','FontSize',14)
fig = gcf;
set(fig,'PaperPositionMode','auto');
fig_pos = get(fig,'PaperPosition');
set(fig,'PaperSize',[fig_pos(3) fig_pos(4)]);
print(fig,'plot_unitaleig_10000','-dpdf')