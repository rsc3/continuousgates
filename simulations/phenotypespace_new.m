function [p,pr] = phenotypespace_new
%calculates the allowed and forbidden gates from a model of gene
%regulation presented in Mayo 2006 (see below)... here we have extended the
%model to allow for two repressors and one activator, and ask which regions
%of parameter space may be accessed which could not be with only a single
%repressor and activator

N = 10000; % number of phentyptes to consider (total)
n = ceil(power(N, 1/4)); %since we have 4 free parameters
randsamp = 1; %pick N points randomly to compare the distributions
muspace = 0;
fullspace = 0;
sepplots = 1; %show separate plots for each combination of randomly sampled plots

%the fixed parameters from A Mayo
alpha = 1.13
gamma = 0.0067
eta = 16.5

%consider parameter ranges over six orders of magnitude
a = logspace(-4,4,n);
c = logspace(-5,5,n);
d = logspace(-4,4,n);
e = logspace(-5,5,n);

%random sets of parameters
if randsamp
    for ind=1:floor(N)
        par = floor(rand(1,4)*n)+1;
        p(:,ind,1) = pi1(alpha, gamma, eta, a(par(1)), c(par(2)), d(par(3)), e(par(4)));
        %now consider the other repressor present (this is not
        %symmetric!)
        p(:,ind,2) = pi2(alpha, gamma, eta, a(par(1)), c(par(2)), d(par(3)), e(par(4)));
        %now consider the co-repression case
        %with the activator present
        p(:,ind,3) = pi3(alpha, gamma, eta, a(par(1)), c(par(2)), d(par(3)), e(par(4)));
        %and without
        p(:,ind,4) = pi4(alpha, gamma, eta, a(par(1)), c(par(2)), d(par(3)), e(par(4)));
        
         %look at three input pispace
        mu(:,ind) = mu_space(alpha, gamma, eta, a(par(1)), c(par(2)), d(par(3)), e(par(4)));
    end
    if sepplots
        subplot(2,2,1);
        plot3(p(3,:,1),p(1,:,1),p(2,:,1),'b.')
        xlabel('\pi_3');
        ylabel('\pi_1');
        zlabel('\pi_2');
        axis([0 1 0 1 0 1]);
        box on
        set(gca,'YDir','reverse');
        title('activation and single repression');
        
        subplot(2,2,2);
        plot3(p(3,:,2),p(1,:,2),p(2,:,2),'g.')
        xlabel('\pi_3');
        ylabel('\pi_1');
        zlabel('\pi_2');
        axis([0 1 0 1 0 1]);
        box on
        set(gca,'YDir','reverse');
        title('activation and repression with a second repressor present');
        
        subplot(2,2,3);
        plot3(p(3,:,3),p(1,:,3),p(2,:,3),'r.')
        xlabel('\pi_3');
        ylabel('\pi_1');
        zlabel('\pi_2');
        axis([0 1 0 1 0 1]);
        box on
        set(gca,'YDir','reverse');
        title('co-repression in the presence of an activator')
        
        subplot(2,2,4);
        plot3(p(3,:,4),p(1,:,4),p(2,:,4),'m.')
        xlabel('\pi_3');
        ylabel('\pi_1');
        zlabel('\pi_2');
        axis([0 1 0 1 0 1]);
        box on
        set(gca,'YDir','reverse');
        title('simpleco-repression');
    end
    
    p2 = sort(p);
    figure;
%     subplot(2,2,1)
%     scatter(p2(3,:,1),p2(2,:,1),sqrt(1./p2(1,:,1))*10,'b.');
%     set(gca,'yscale','log','xscale','log');
%     title('activation and single repression');
%     subplot(2,2,2)
    scatter(p2(3,:,2),p2(2,:,2),sqrt(1./p2(1,:,2))*10,'g.');
%     set(gca,'yscale','log','xscale','log');
%     title('activation and repression with a second repressor present');
%     subplot(2,2,3)
hold on
    scatter(p2(3,:,3),p2(2,:,3),sqrt(1./p2(1,:,3))*10,'r.');
    hold off
    set(gca,'yscale','log','xscale','log');
%     title('co-repression in the presence of an activator')
%     subplot(2,2,4)
%     scatter(p2(3,:,4),p2(2,:,4),1./p2(1,:,4),'m.')
%     set(gca,'yscale','log','xscale','log');
%     title('simple co-repression');

    figure;
    plot3(p(3,:,1),p(1,:,1),p(2,:,1),'b.',p(3,:,2),p(1,:,2),p(2,:,2),'g.',p(3,:,3),p(1,:,3),p(2,:,3),'r.',p(3,:,4),p(1,:,4),p(2,:,4),'m.','MarkerSize',1)
    legend({'A + R1, R2 absent','A + R1, R2 present','R1 + R2, A present','R1 + R2, A absent'});
    xlabel('\pi_3');
    ylabel('\pi_1');
    zlabel('\pi_2');
    axis([0 1 0 1 0 1]);
    box on
    set(gca,'YDir','reverse');
end

if muspace
    figure;
    scatter3(mu(2,:),mu(3,:),mu(4,:),1./mu(1,:))
    xlabel('\mu_2');
    ylabel('\mu_3');
    zlabel('\mu_4');
    axis([0 1 0 1 0 1]);
    box on
    figure
    scatter(mu(2,:),mu(3,:),1./mu(1,:))
end


%consider all combinations of parameters

if fullspace
    ind = 0;
    for i=1:n
        for j=1:n
            for k=1:n
                for l=1:n
                    ind = ind+1;
                    %first replicate the phenotypes from Mayo et al
                    p(:,ind,1) = pi1(alpha, gamma, eta, a(i), c(j), d(k), e(l));
                    %now consider the other repressor present (this is not
                    %symmetric!)
                    p(:,ind,2) = pi2(alpha, gamma, eta, a(i), c(j), d(k), e(l));
                    %now consider the co-repression case
                    %with the activator present
                    p(:,ind,3) = pi3(alpha, gamma, eta, a(i), c(j), d(k), e(l));
                    %and without
                    p(:,ind,4) = pi4(alpha, gamma, eta, a(i), c(j), d(k), e(l));
                end
            end
        end
    end
    figure;
    plot3(p(3,:,1),p(1,:,1),p(2,:,1),'k.',p(3,:,2),p(1,:,2),p(2,:,2),'b.',p(3,:,3),p(1,:,3),p(2,:,3),'g.',p(3,:,4),p(1,:,4),p(2,:,4),'r.')
    xlabel('\pi_3');
    ylabel('\pi_1');
    zlabel('\pi_2');
    legend({'A + R1, R2 absent','A + R1, R2 present','R1 + R2, A present','R1 + R2, A absent'});
    set(gca,'YDir','reverse');
    axis([0 1 0 1 0 1]);
    box on
    
    %consider the plateau logic gates
    pg = [1/100 1/100 10/100; 1/100 10/100 10/100; 1/100 10/100 100/100; 1/1000 10/1000 100/1000];
    hold on
    plot3(pg(:,1),pg(:,2),pg(:,3),'*','MarkerSize',10)
    hold off
end

function activity = P(alpha, gamma, eta, A, R, S, a, c, d, e)
%this is the CRIF of three inputs
%A, R, S are scalar input variables
%alpha, gamma, eta are scalar internal parameters
%a, c, d, e are (vector) external parameters
activity = (alpha - gamma) * a .* (1 + eta * A * d) ./ (1 + a + (a*eta + 1) .* d * A + c * R + e * S) + gamma;

%these functions define the 8 different plateaus achieved from the three
%inputs
function plateau = P1(alpha, gamma, eta, a, c, d, e)
plateau = P(alpha, gamma, eta, 0, 1, 0, a, c, d, e);
function plateau = P2(alpha, gamma, eta, a, c, d, e)
plateau = P(alpha, gamma, eta, 0, 0, 0, a, c, d, e);
function plateau = P3(alpha, gamma, eta, a, c, d, e)
plateau = P(alpha, gamma, eta, 1, 1, 0, a, c, d, e);
function plateau = P4(alpha, gamma, eta, a, c, d, e)
plateau = P(alpha, gamma, eta, 1, 0, 0, a, c, d, e);
function plateau = P5(alpha, gamma, eta, a, c, d, e)
plateau = P(alpha, gamma, eta, 0, 1, 1, a, c, d, e);
function plateau = P6(alpha, gamma, eta, a, c, d, e)
plateau = P(alpha, gamma, eta, 0, 0, 1, a, c, d, e);
function plateau = P7(alpha, gamma, eta, a, c, d, e)
plateau = P(alpha, gamma, eta, 1, 1, 1, a, c, d, e);
function plateau = P8(alpha, gamma, eta, a, c, d, e)
plateau = P(alpha, gamma, eta, 1, 0, 1, a, c, d, e);

%the phenotypes due to a single repressor and activator with the second
%repressor absent
function ratios = pi1(alpha, gamma, eta, a, c, d, e)
ratios = [P1(alpha, gamma, eta, a, c, d, e)./P4(alpha, gamma, eta, a, c, d, e);
    P2(alpha, gamma, eta, a, c, d, e)./P4(alpha, gamma, eta, a, c, d, e);
    P3(alpha, gamma, eta, a, c, d, e)./P4(alpha, gamma, eta, a, c, d, e)];

%the phenotypes due to a single repressor and activator with the second
%repressor present
function ratios = pi2(alpha, gamma, eta, a, c, d, e)
ratios = [P5(alpha, gamma, eta, a, c, d, e)./P8(alpha, gamma, eta, a, c, d, e);
    P6(alpha, gamma, eta, a, c, d, e)./P8(alpha, gamma, eta, a, c, d, e);
    P7(alpha, gamma, eta, a, c, d, e)./P8(alpha, gamma, eta, a, c, d, e)];

%the phenotypes due to co-repression in the presence of the activator
function ratios = pi3(alpha, gamma, eta, a, c, d, e)
ratios = [P7(alpha, gamma, eta, a, c, d, e)./P4(alpha, gamma, eta, a, c, d, e);
    P3(alpha, gamma, eta, a, c, d, e)./P4(alpha, gamma, eta, a, c, d, e);
    P8(alpha, gamma, eta, a, c, d, e)./P4(alpha, gamma, eta, a, c, d, e)];

%the phenotypes due to co-repression in the absence of the activator
function ratios = pi4(alpha, gamma, eta, a, c, d, e)
ratios = [P5(alpha, gamma, eta, a, c, d, e)./P2(alpha, gamma, eta, a, c, d, e);
    P1(alpha, gamma, eta, a, c, d, e)./P2(alpha, gamma, eta, a, c, d, e);
    P6(alpha, gamma, eta, a, c, d, e)./P2(alpha, gamma, eta, a, c, d, e)];

%calculate the 3-input naive phenotypes

function tp = mu_space(alpha, gamma, eta, a, c, d, e)
    %the range
    tp = [P(alpha, gamma, eta, 0, 1, 1, a, c, d, e) ./ P(alpha, gamma, eta, 1, 0, 0, a, c, d, e);
        %the activator
        P(alpha, gamma, eta, 1, 1, 1, a, c, d, e) ./ P(alpha, gamma, eta, 1, 0, 0, a, c, d, e);
        %repressor one
        P(alpha, gamma, eta, 0, 0, 1, a, c, d, e) ./ P(alpha, gamma, eta, 1, 0, 0, a, c, d, e);
        %repressor two
        P(alpha, gamma, eta, 0, 1, 0, a, c, d, e) ./ P(alpha, gamma, eta, 1, 0, 0, a, c, d, e);];