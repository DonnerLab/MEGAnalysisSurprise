function plot_binned_scatter(x,y,nbins,fittype,c)

% x = vector of x inputs
% y = vector of x inputs
% nbins = desired number of bins (sorted by x)
% fittype = regression fit type; 'linear' or 'logistic'
% c = colour of points in scatterplot (in rgb)

breaks = 0:length(x)/nbins:length(x);
breaks(1) = 1; if breaks(end)~=length(x), breaks(end+1)=length(x); end

[~,sort_x] = sort(x);
x_b=[]; y_b=[]; y_SE=[];
for bin = 1:length(breaks)-1;
    x_b(length(x_b)+1,1) = mean(x(sort_x(ceil(breaks(bin)):floor(breaks(bin+1)))));
    y_b(length(y_b)+1,1) = mean(y(sort_x(ceil(breaks(bin)):floor(breaks(bin+1)))));
    y_SE(length(y_SE)+1,1) = std(y(sort_x(ceil(breaks(bin)):floor(breaks(bin+1)))))./sqrt(length(y(sort_x(ceil(breaks(bin)):floor(breaks(bin+1))))));
end

yrange=[min(y_b-y_SE)-(min(abs(diff(y_b)))) max(y_b+y_SE)+(min(abs(diff(y_b))))];
xrange=[min(x_b)-(min(abs(diff(x_b)))) max(x_b)+(min(abs(diff(x_b))))];

plot([0 0],[yrange],'LineStyle','--','Color',[0.75 0.75 0.75])
plot([xrange],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75])

if strcmp(fittype,'linear')
    m = regstats(y,x,'linear',{'beta'});
    plot([min(x_b) max(x_b)],[m.beta(1)+(m.beta(2)*min(x_b)) m.beta(1)+(m.beta(2)*max(x_b))],'k')
elseif strcmp(fittype,'logistic')
    m = glmfit(x, [y ones(size(y))], 'binomial');
    yvals = [1/(1+exp(-(m(1)+(m(2)*min(x_b))))) 1/(1+exp(-(m(1)+(m(2)*max(x_b)))))];
    plot([min(x_b) max(x_b)],yvals,'k')
end

for b=1:length(x_b), plot([x_b(b) x_b(b)],[y_b(b)-y_SE(b) y_b(b)+y_SE(b)],'LineWidth',0.5,'Color',[0.65 0.65 0.65]), end
s = scatter(x_b,y_b,45,c,'o','filled'); set(s,'MarkerFaceColor',c,'MarkerEdgeColor',[0 0 0]);
ylim(yrange), xlim(xrange)

