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
A(:,m) = mean(A(:,1:spark-1),2);% Introducing spark though unneccessary
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
x_k1= [];
t   = 0.5;
tic
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
toc
%% Results
x_wmp = zeros(m,1);
x_wmp(S_k) =   x_k;
disp(['No of non-zeros in original solution: ' num2str(K)]);
disp(['OMP provided a solution with : ' num2str(length(S_k)) ' non zero elements']);
disp(['||r_k||  : ' num2str(norm(r_k)) ' and error : ' num2str(norm(x-x_wmp))]);
