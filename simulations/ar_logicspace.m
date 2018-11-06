function p = phenotypespace

goodRR = 0;

N = 1000; % number of points to consider (total)
B = 36; %number of extreme gates to analyze
n = ceil(power(N, 1/3));
e = ceil(power(B, 1/2));
plotit = 1;






%% %%%%%%%%%%%MAYO 2006 MODEL PARAMETERS%%%%%%%%%%%%%%%
alpha = 1.13 * ones(1,n)
gamma = 0.0067 * ones(1,n)
eta = 16.5 * ones(1,n)

a = logspace(-3,3,n);
c = logspace(-3,3,n);
d = logspace(-3,3,n);

Kr = 1.3; %uM
Ka = 2;%mM
nr = 4;
na = 2;


%% %%%%%%%%SETTY AND OR GATES%%%%%%%%%%%%%%%%%%%
Kr = 1.2;
Ka = 1.8;
Vand = [1 70 2000 17 7000];
Vor  = [10 1700 15  400 50];

    Ai=logspace(-2,4);
    Ri=logspace(-1,2);

    li = length(Ai);
    A = fliplr(hill01(Ai,Ka,na));
    R = hill01(Ri,Kr,nr);
    AM=repmat(A,length(R),1);
    RM=repmat(R,length(A),1);
    
AND=f(Vand,AM,RM');
OR =f(Vor,AM,RM');

[ra,aa,la]=logicspace([AND(1,1) AND(1,end) AND(end,1) AND(end,end)]);
[ro,ao,lo]=logicspace([OR(1,1) OR(end,1) OR(1,end) OR(end,end)]);
%%% %%%%%%%%%%%%NEW PARAMETERS
a = logspace(-2,1,n);
c = logspace(2,4,n);
d = logspace(-1,2,n);
gamma(:)=0

ind = 0;
for i=1:n
    for j=1:n
        for k=1:n
            ind = ind+1;
            p(ind,:) = pi(alpha(i), gamma(j), eta(k), a(i), c(j), d(k));
            param(ind,:) = [a(i) c(j) d(k)];
            %the distance to the perfect binary gates of phentype space
            b(ind,:) = sort([P1(alpha(i), gamma(j), eta(k), a(i), c(j), d(k))
                P2(alpha(i), gamma(j), eta(k), a(i), c(j), d(k))
                P3(alpha(i), gamma(j), eta(k), a(i), c(j), d(k))
                P4(alpha(i), gamma(j), eta(k), a(i), c(j), d(k))]);
        end
    end
end

[s,I]=sort(p(:,1));
extgateparams = param(I(1:B),:);
[r,a,l]=logicspace(b);

if plotit
    figure(104);clf;
    hl1 = line([0 0.5 1 0],[0 1 0 0]);
    set(hl1,'color','k','linestyle',':');
    ylim([0 1]);
    hold on;
    scatter(l,a,r*100,'MarkerEdgeColor','k')
    xlabel('\leftarrow more OR   more AND \rightarrow','FontSize',18);
    ylabel('asymmetry', 'FontSize',18);
    axis0;
    set(gca,'xtick',[0 .25 .5 .75 1],'FontSize',18);
    hold off;
end

%find the gates with high RR
if goodRR
    goodRR = find(p(:,1) < .01);
    figure;
    plot3(param(:,1),param(:,2),param(:,3),'.');
    hold on
    plot3(param(goodRR,1),param(goodRR,2),param(goodRR,3),'r*')
    hold off
end
%find the most extreme gates

%[s, I] = sort(dis);
%extgateparams = param(I(1:B),:);
figure
for i=1:B
    subplot(e,e,i);
    Ai=logspace(-2,4);
    Ri=logspace(-1,2);
    A = fliplr(hill01(Ai,Ka,na));
    R = hill01(Ri,Kr,nr);
    ae = extgateparams(i,1);
    ce = extgateparams(i,2);
    de = extgateparams(i,3);
    for j=1:length(A)
        for k=1:length(R)
            CRIF(j,k) = P(alpha(1), gamma(1), eta(1), A(j), R(k),ae,ce,de);
        end
    end
    surf(CRIF,'EdgeColor','none');
    set(gca,'xtick',0,'ytick',0,'ztick',0)
    rotate3d on;
end

figure;
scatter3(param(1:10:end,1),param(1:10:end,2),param(1:10:end,3),'k.');
hold on;
scatter3(extgateparams(:,1),extgateparams(:,2),extgateparams(:,3),'ro','filled');
hold off;
xlabel('a')
ylabel('c')
zlabel('d')
set(gca,'xscale','log','yscale','log','zscale','log')

function activity = P(alpha, gamma, eta, A, R, a, c, d)
%alpha, gamma, eta are scalar internal parameters
%A, R are scalar input variables
%a, c, d are vector external parameters
activity = (alpha - gamma) * a .* (1 + eta * A * d) ./ (1 + a + (a*eta + 1) .* d * A + c * R) + gamma;


function activity = f(V,A,R)
%same model, different parameters V1-V5
activity = V(1)*(1+V(2)*A+V(3)*R)./(1+V(4)*A+V(5)*R)

function plateau = P1(alpha, gamma, eta, a, c, d)
plateau = P(alpha, gamma, eta, 0, 1, a, c, d);

function plateau = P2(alpha, gamma, eta, a, c, d)
plateau = P(alpha, gamma, eta, 0, 0, a, c, d);

function plateau = P3(alpha, gamma, eta, a, c, d)
plateau = P(alpha, gamma, eta, 1, 1, a, c, d);

function plateau = P4(alpha, gamma, eta, a, c, d)
plateau = P(alpha, gamma, eta, 1, 0, a, c, d);

function [r,a,l] = logicspace(b)
%transforms binary gate data into logic space data

r = log10(b(:,4)./b(:,1));
beta = log10(b(:,4)./b(:,2)) ./ r;
alpha = log10(b(:,4)./b(:,3)) ./ r;

l = (alpha + beta)/2;
a = (beta - alpha);

function ratios = pi(alpha, gamma, eta, a, c, d)
ratios = [P1(alpha, gamma, eta, a, c, d)./P4(alpha, gamma, eta, a, c, d);
    P2(alpha, gamma, eta, a, c, d)./P4(alpha, gamma, eta, a, c, d);
    P3(alpha, gamma, eta, a, c, d)./P4(alpha, gamma, eta, a, c, d)];

function p = hill01(x,affinity, cooperativity)
p = 1./(1 + (x/affinity).^cooperativity);