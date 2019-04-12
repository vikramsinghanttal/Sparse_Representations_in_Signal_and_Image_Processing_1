% ref : http://math.mit.edu/~liewang/OMP.pdf
% Orthogonal Matching Pursuit for Sparse Signal Recovery With Noise
%% Complexity Analysis
% Most complex operations are :
% A'*r_k       ~ O(mn) 
% inv(As'*As)  ~ O(K*K*m)
% Overall complexity of OMP is O(mnk)


clear all;
clc;
% rng(120);

%% Setting the parameters
n = 50; m =120;
threshold = 0.1;
montecarloiterations = 1000;
error_omp = zeros(15,1);
Pe_sup_omp= zeros(15,1);
for mci=1:montecarloiterations
    for K=1:15

    %% Generating A, x and b
    % K = randi([1,n/2]);             % No. of nonzero parameters in x
    nonz_idx = randi([1,120],K,1);  % Indices which will contain the non zero elements in x

    % A=sqrt(0.5)*(randn(n,m)+1i*randn(n,m));
    A = randn(n,m);                 % Dictionary matrix
    % A  = A./vecnorm(A);           % Not supported by R2015a
    A = A*diag(1./sqrt(diag(A'*A)));% making columns unit norm
    x = zeros(m,1);
    x(nonz_idx) = 1+rand(K,1); 
    b = A*x;
        %% Auxilliary Variables
        k   = 0;
        r_k = b;
        S_k = [];
        E_i = inf;
        z_i = 0;
        [A_sorted,sorted_idx] = sort((A'*b));
%         tic
        %% OMP Algorithm
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
    end
    mci
end
figure(1)
plot(error_omp/montecarloiterations);
title('L-2 Error vs K');
ylabel('L-2 error in estimate');
xlabel('K');
p.LineWidth = 2;
lgd1.FontSize = 14;

figure(2)
plot(Pe_sup_omp/montecarloiterations);
title('P_e Support vs K');
ylabel('Prob of mismatch in support');
xlabel('K');
p.LineWidth = 2;
lgd1.FontSize = 14;