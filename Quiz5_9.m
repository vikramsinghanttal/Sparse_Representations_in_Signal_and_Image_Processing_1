% Quiz 5: Greedy and Relexation Methods 
clear all;
clc;
A = [ 0.1817   0.5394 -0.1197  0.6404;
      0.6198   0.1994  0.0946 -0.3121;
     -0.7634  -0.8181  0.9883  0.7018];
% A = A./vecnorm(A);
b = [1.1862; -0.1158; -0.1093];
threshold = 0.1;
m = size(A,2);
n = size(A,2);
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
while(norm(r_k)>threshold && k<2)
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
    fprintf(['[WMP] Residual ' num2str(norm(r_k)) ' at for column ' num2str(i) '\n']);
end