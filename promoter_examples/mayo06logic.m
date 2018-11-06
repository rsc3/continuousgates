
pi = [.07 .01 .003 .007 .01 .01 .01 .07 .01 .01 .06 .05 .05; 
    .14 .15 .08 .20 .14 .053 .092 .32 .009 .009 .23 .3 .11;
    .28 .26 .11 .15 .25 .4 .09 .63 .28 .11 .41 .18 .081];

pi = sort(pi);


r = 1./pi(1,:);

a = log10(pi(3,:)./pi(2,:))./log10(r);

l = -log10(pi(2,:).*pi(3,:))./(2*log10(r));

%log10(sqrt((1./pi(2,:)).*(1./pi(3,:))))./log10(r);

scatter(l,a,log10(r).^2*100,'MarkerEdgeColor','k','linewidth',2)
hl1 = line([0 0.5 1 0],[0 1 0 0]);
set(hl1,'color','k','linestyle',':');
hold on;

scatter([0.15 0.15 0.15],[.5 .7 .9],[1 4 9]*100,'linewidth',2)
text([0.1 0.1 0.1]-.05,[.5 .7 .9],{'10^1','10^2','10^3'},'FontSize',18)
    
%show the wild type in red
scatter(l(1),a(1),log10(r(1)).^2*100,'MarkerEdgeColor','r','linewidth',2)

%show the CRP dominant clones in blue
scatter(l([4 12 13]),a([4 12 13]),log10(r([4 12 13])).^2*100,'MarkerEdgeColor','b','linewidth',2)

%plot dashed boundaries
hbound1 = line([0.125 0.875],[0.25 0.25])
set(hbound1,'color','k','linestyle','--')
hbound2 = line([0.25 0.25],[0 0.25])
set(hbound2,'color','k','linestyle','--')
hbound3 = line([0.75 0.75],[0 0.25])
set(hbound3,'color','k','linestyle','--')

hbound4 = line([0.375 0.625],[0.75 0.75])
set(hbound4,'color','k','linestyle','--')
hbound5 = line([0.375 0.375],[0.25 0.75])
set(hbound5,'color','k','linestyle','--')
hbound6 = line([0.625 0.625],[0.25 0.75])
set(hbound6,'color','k','linestyle','--')

hold off;
set(gca,'xtick',[0 .25 .5 .75 1],'ytick',[0 .25 .5 .75 1],'FontSize',18,'xlim',[0 1],'ylim',[0 1]);
xlabel('\leftarrow more OR   more AND \rightarrow \it{(l)}');
ylabel('asymmetry \it{(a)}');
axis0;