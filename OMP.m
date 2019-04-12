clear all;
clc
% rng(120);
n=50; m=150;
spark = 13;
% A=sqrt(0.5)*(randn(n,m)+1i*randn(n,m));
A = randn(n,m);
A(:,m) = mean(A(:,1:spark-1),2);

%% Normalizing the columns of A
A  = A./vecnorm(A);
threshold = 0.1;
K = randi([1,n/2]);
nonz_idx = randi([1,120],K,1);

A = randn(n,m);
b = randn(n,1);
x = zeros(m,1);
x(nonz_idx) = randn(K,1); 

k=0;
r_k = b;
S_k = [];

E_i = inf;
z_i = 0;

while(norm(r_k)>threshold && k<m)
    z   = A'*r_k;
    [M,I]= max(abs(z));
    S_k = [S_k I];
    As  = A(:,S_k);
    x_k = (As'*As)\(As'*b);
    r_k = b - As*x_k;
    k=k+1;
end

x_est = zeros(m,1);
x_est(S_k)=x_k;
disp(['No of non-zeros in original solution: ' num2str(K)]);
disp(['OMP provided a solution with : ' num2str(length(S_k)) ' non zero elements']);
disp(['||r_k||  : ' num2str(norm(r_k)) ' and error : ' num2str(norm(x-x_est))]);