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
montecarloiterations = 10;
error_lsomp = zeros(15,1);
Pe_sup_lsomp= zeros(15,1);
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
        A_idx = 1:1:m; 
        best_idx = 0;

    %% LS-OMP Algorithm
        while(E_i > threshold && k<m)
            E_i = inf;
            for i  = 1:m-k
                As = [A(:,S_k) A(:,A_idx(i))];
                x_i= (As'*As)\(As'*b);
                norm_rk=norm(As*x_i-b);
                if E_i > norm_rk
                   E_i = norm_rk;
                   best_idx  = i;
                   x_k = x_i;
                end
            end
            k = k+1;
            S_k = [S_k A_idx(best_idx)];
            A_idx(best_idx)=[];
        end

    %% Results
        x_lsomp = zeros(m,1);
        x_lsomp(S_k)=x_k;
        error_lsomp(K) = error_lsomp(K)+norm(x_lsomp-x)/norm(x);
        Pe_sup_lsomp(K)= Pe_sup_lsomp(K) + 1-sum(x&x_lsomp)/max(nnz(x),nnz(x_lsomp));
    end
    mci
end
figure(1)
plot(error_lsomp/montecarloiterations);
title('L-2 Error vs K');
ylabel('L-2 error in estimate');
xlabel('K');
p.LineWidth = 2;
lgd1.FontSize = 14;

figure(2)
plot(Pe_sup_lsomp/montecarloiterations);
title('P_e Support vs K');
ylabel('Prob of mismatch in support');
xlabel('K');
p.LineWidth = 2;
lgd1.FontSize = 14;