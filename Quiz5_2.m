% Quiz 5: Greedy and Relexation Methods 
clear all;
clc;
A = [ 0.1817   0.5394 -0.1197  0.6404;
      0.6198   0.1994  0.0946 -0.3121;
     -0.7634  -0.8181  0.9883  0.7018];
% A = A./vecnorm(A);
b = [1.1862; -0.1158; -0.1093];
S_k = [];
r_k = b;

%% First Iteration
[M,I] = max(abs(A'*r_k));
fprintf(['[OMP] max correlation ' num2str(M) ' at for column ' num2str(I) '\n']);
S_k = [S_k, I];
As = A(:,S_k);
x_k = pinv(As)*b;
r_k = b - As*x_k;
fprintf(['[OMP] Residual ' num2str(norm(r_k)) '\n']);

%% Second Iteration
[M,I] = max(abs(A'*r_k));
fprintf(['[OMP] max correlation ' num2str(M) ' at for column ' num2str(I) '\n']);
S_k = [S_k, I];
As = A(:,S_k);
x_k = pinv(As)*b;
r_k = b - As*x_k;
fprintf(['[OMP] Residual ' num2str(norm(r_k)) '\n']);