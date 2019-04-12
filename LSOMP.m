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
spark = 13;
threshold = 0.1;

%% Generating A, x and b
K = randi([1,n/2]);             % No. of nonzero parameters in x
nonz_idx = randi([1,120],K,1);  % Indices which will contain the non zero elements in x

% A=sqrt(0.5)*(randn(n,m)+1i*randn(n,m));
A = randn(n,m);                 % Dictionary matrix
% A(:,m) = mean(A(:,1:spark-1),2);% Introducing spark though unneccessary
% A  = A./vecnorm(A);           % Not supported by R2015a
A = A*diag(1./sqrt(diag(A'*A)));% making columns unit norm
b = randn(n,1);
x = zeros(m,1);
x(nonz_idx) = randn(K,1); 

%% Auxilliary Variables
k   = 0;
r_k = b;
S_k = [];
E_i = inf;
z_i = 0;
A_idx = 1:1:m; 
best_idx = 0;
tic
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
toc
%% Results
x_lsomp = zeros(m,1);
x_lsomp(S_k)=x_k;
r_k=norm(A*x_lsomp-b);
disp(['No of non-zeros in original solution: ' num2str(K)]);
disp(['OMP provided a solution with : ' num2str(length(S_k)) ' non zero elements']);
disp(['||r_k||  : ' num2str((r_k)) ' and error : ' num2str(norm(x-x_lsomp))]);
