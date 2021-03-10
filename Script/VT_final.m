close all
hfig = figure(3);
rng default;
x = rand([1 400]);
y = rand([1 400]);
voronoi(x,y)
xlim([0 1])
ylim([0 1])
% axis equal
box on
ax = gca;
ax.XColor = 'black';
ax.LineWidth = 1.5
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
set(gca,'linewidth', 2);
% figWidth = 64;
% figHight = 64;
% set(hfig,'PaperUnits','centimeters');
% set(hfig,'PaperPosition',[0 0 figWidth figHight])
% fileout = [mat2str(1)];
% print(hfig,[fileout,'fcc'],'-r300','-dpng')
