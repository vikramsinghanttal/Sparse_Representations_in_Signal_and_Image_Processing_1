clear all;
clc;
% rng(120);

%% Setting the parameters
n = 50; m =120;
threshold = 0.1;
montecarloiterations = 10;
x_lower = 1;
x_upper = 1;
error_IRLS = zeros(15,1);
Pe_sup_IRLS= zeros(15,1);
error_lsomp = zeros(15,1);
Pe_sup_lsomp= zeros(15,1);
error_mp = zeros(15,1);
Pe_sup_mp= zeros(15,1);
error_wmp = zeros(15,1);
Pe_sup_wmp= zeros(15,1);
error_thr = zeros(15,1);
Pe_sup_thr= zeros(15,1);
for mci=1:montecarloiterations
    for K=1:15
        %% Generating A, x and b
        nonz_idx = randi([1,120],K,1);  % Indices which will contain the non zero elements in x
        
        A = randn(n,m);                 % Dictionary matrix
        % A  = A./vecnorm(A);           % Not supported by R2015a
        A = A*diag(1./sqrt(diag(A'*A)));% making columns unit norm
        x = zeros(m,1);
        x(nonz_idx) = x_lower+(x_upper-x_lower)*rand(K,1); 
        b = A*x;
        %% Auxilliary Variables
        k   = 0;
        r_k = b;
        S_k = [];
        E_i = inf;
        z_i = 0;
        x_k = ones(m,1);

        %% OMP Algorithm
        while(norm(r_k)>threshold && k<5)
            rho_k = sqrt(abs(x_k));
            w_k = rho_k./(x_k.^2);
            W_k = diag(1./w_k);
            WA_k = W_k*A';
            x_k = (WA_k/(A*WA_k))*b;
            k = k+1;
            r_k = b - A*x_k;
        end

            %% Results
        x_IRLS = x_k;
%         x_IRLS(S_k)=x_k;
        error_IRLS(K) = error_IRLS(K)+norm(x_IRLS-x)/norm(x);
        Pe_sup_IRLS(K)=Pe_sup_IRLS(K)+(1-sum(x&x_IRLS)/max(nnz(x),nnz(x_IRLS)));

    end
    mci
end
error_IRLS   = error_IRLS/montecarloiterations;
% error_lsomp = error_lsomp/montecarloiterations;

Pe_sup_IRLS  = Pe_sup_IRLS/montecarloiterations;
% Pe_sup_lsomp= Pe_sup_lsomp/montecarloiterations;
x_ax = (1:1:15)';
p = figure(1); clf;
set(p,'Position',[115 100 600 400]);
% set(gca,'xscale','log')
hold on
p = semilogx(x_ax,error_IRLS,'r');
q.LineWidth = 1.5;
title(['L_2 Estimation Error vs K (Error Tolerance = ' num2str(threshold) ' )']);
ylabel('Avg. Relative L_2 error in estimates');
xlabel('Cordinality of the solution (K)');

q = figure(2); clf;
set(q,'Position',[730 100 600 400]);
hold on
q = semilogx(x_ax,Pe_sup_IRLS,'r');
q.LineWidth = 1.5;
title(['P_e Support vs K (Error Tolerance = ' num2str(threshold) ' )']);
ylabel('Probability of mismatch in support');
xlabel('Cordinality of the Solution (K)');