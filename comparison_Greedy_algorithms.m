clear all;
clc;
% rng(120);

%% Setting the parameters
n = 50; m =120;
threshold = 0.1;
error_omp = zeros(15,1);
Pe_sup_omp= zeros(15,1);

error_lsomp = zeros(15,1);
Pe_sup_lsomp= zeros(15,1);

error_mp = zeros(15,1);
Pe_sup_mp= zeros(15,1);

error_wmp = zeros(15,1);
Pe_sup_wmp= zeros(15,1);

error_th = zeros(15,1);
Pe_sup_th= zeros(15,1);
montecarloiterations = 10;
for mci = 1:montecarloiterations
    for K=1:15                          % No. of nonzero parameters in x
        %% Generating A, x and b
        nonz_idx = randi([1,120],K,1);  % Indices which will contain the non zero elements in x

        % A=sqrt(0.5)*(randn(n,m)+1i*randn(n,m));
        A = randn(n,m);                 % Dictionary matrix
        % A  = A./vecnorm(A);           % Not supported by R2015a
        A = A*diag(1./sqrt(diag(A'*A)));% making columns unit norm
        x = zeros(m,1);
        x(nonz_idx) = randn(K,1);
        b = A*x;

        %% Auxilliary Variables
        k   = 0;
        r_k = b;
        S_k = [];
        E_i = inf;
        z_i = 0;


        %% OMP Algorithm
%         tic
        while(norm(r_k)>threshold && k<m)
            z   = A'*r_k;
            [M,I]= max(abs(z));
            S_k = [S_k I];
            As  = A(:,S_k);
            x_k = (As'*As)\(As'*b);
            r_k = b - As*x_k;
            k=k+1;
        end
%         toc
    %% Results
        x_omp = zeros(m,1);
        x_omp(S_k)=x_k;
        error_omp(K) = error_omp(K)+norm(x_omp-x)/norm(x);
        Pe_sup_omp(K)=Pe_sup_omp(K)+(1-sum(x&x_omp)/max(nnz(x),nnz(x_omp)));


        %% Auxilliary Variables
        k   = 0;
        r_k = b;
        S_k = [];
        E_i = inf;
        z_i = 0;
        A_idx = 1:1:m; 
        best_idx = 0;
%         tic
        %% LS-OMP Algorithm
        while(E_i > threshold && k<m)
            E_i = inf;
            for i  = 1:m-k
                As = [A(:,S_k) A(:,A_idx(i))];
                x_i= (As'*As)\As'*b;
                r_k=norm(As*x_i-b);
                if E_i > r_k
                   E_i = r_k;
                   ii  = i;
                   x_k = x_i;
                   best_idx = A_idx(i);
                end
            end
            k = k+1;
            S_k = [S_k best_idx];
            A_idx(ii)=[];
        end
%         toc
        %% Results
        x_lsomp = zeros(m,1);
        x_lsomp(S_k)=x_k;
        error_lsomp(K) = error_lsomp(K)+norm(x_lsomp-x)/norm(x);
        Pe_sup_lsomp(K)= Pe_sup_lsomp(K) + 1-sum(x&x_lsomp)/max(nnz(x),nnz(x_lsomp));


        %% Auxilliary Variables
        k   = 0;
        r_k = b;
        S_k = [];
        E_i = inf;
        z_i = 0;
        x_k1 =[];
%         tic
        %% Matching Pursuit Algorithm
        while(norm(r_k)>threshold && k<m)
            z   = A'*r_k;
            [M,I]= max(abs(z));
            S_k = [S_k I];
            As  = A(:,S_k);
            x_k = [x_k1;A(:,I)'*r_k];
            r_k = b - As*x_k;
            k   = k+1;
            x_k1= x_k;
        end
%         toc
        %% Results
        x_mp = zeros(m,1);
        x_mp(S_k)=x_k;

        error_mp(K) = error_mp(K)+norm(x_mp-x)/norm(x);
        Pe_sup_mp(K)= Pe_sup_mp(K)+ 1-sum(x&x_mp)/max(nnz(x),nnz(x_mp));

        k   = 0;
        r_k = b;
        S_k = [];
        E_i = inf;
        z_i = 0;
        x_k1= [];
        t   = 0.5;
%         tic
        %% Weak Matching Pursuit Algorithm
        while(norm(r_k)>threshold && k<m)
            for i=1:m
                z_i   = A(:,i)'*r_k;
                if abs(z_i) > t*norm(r_k)
                    break;
                end
            end
            S_k = [S_k i];
            As  = A(:,S_k);
            x_k = [x_k1;A(:,i)'*r_k];
            r_k = b - As*x_k;
            k   = k+1;
            x_k1= x_k;
        end
%         toc
        %% Results
        x_wmp = zeros(m,1);
        x_wmp(S_k) =   x_k;
        error_wmp(K) = error_wmp(K) + norm(x_wmp-x)/norm(x);
        Pe_sup_wmp(K)= Pe_sup_wmp(K)+ 1-sum(x&x_wmp)/max(nnz(x),nnz(x_wmp));

%         disp(['No of non-zeros in original solution: ' num2str(K)]);
%         disp(['OMP provided a solution with : ' num2str(length(S_k)) ' non zero elements']);
%         disp(['||r_k||  : ' num2str(norm(r_k)) ' and error : ' num2str(norm(x-x_omp)/norm(x))]);

    end
    mci
end
figure(1)
hold on
plot(error_omp/montecarloiterations);
plot(error_lsomp/montecarloiterations);
plot(error_mp/montecarloiterations);
plot(error_wmp/montecarloiterations);
title('L-2 Error vs K');
xlabel('L-2 error in estimate');
ylabel('K');
legend('OMP','LS-OMP','MP','WMP');

figure(2)
hold on
plot(Pe_sup_omp/montecarloiterations);
plot(Pe_sup_lsomp/montecarloiterations);
plot(Pe_sup_mp/montecarloiterations);
plot(Pe_sup_wmp/montecarloiterations);
title('P_e Support vs K');
xlabel('Prob of mismatch in support');
ylabel('K');
legend('OMP','LS-OMP','MP','WMP');