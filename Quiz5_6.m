% Quiz 5: Greedy and Relexation Methods 
clear all;
clc;
A = [ 0.1817   0.5394 -0.1197  0.6404;
      0.6198   0.1994  0.0946 -0.3121;
     -0.7634  -0.8181  0.9883  0.7018];
% A = A./vecnorm(A);
b = [1.1862; -0.1158; -0.1093];
t =0.5;

m = size(A,2);
k   = 0;
r_k = b;
S_k = [];
E_i = inf;
z_i = 0;
A_idx = 1:1:m; 
best_idx = 0;

%% First Iteration
E_i = inf;
for i  = 1:m-k
	As = [A(:,S_k) A(:,A_idx(i))];
	x_i= (As'*As)\(As'*b);
    r_k= As*x_i-b;
	norm_rk=norm(r_k);
	if E_i > norm_rk
        E_i = norm_rk;
        best_idx  = i;
        x_k = x_i;
	end
end
fprintf(['[LS-OMP] max correlation ' num2str(E_i) ' at for column ' num2str(A_idx(best_idx)) '\n']);
k = k+1;
S_k = [S_k A_idx(best_idx)];
A_idx(best_idx)=[];

%% Second Iteration
E_i = inf;
for i  = 1:m-k
	As = [A(:,S_k) A(:,A_idx(i))];
	x_i= (As'*As)\(As'*b);
    r_k= As*x_i-b;
	norm_rk=norm(r_k);
	if E_i > norm_rk
        E_i = norm_rk;
        best_idx  = i;
        x_k = x_i;
	end
end
fprintf(['[LS-OMP] Residual ' num2str(E_i) ' at for column ' num2str(A_idx(best_idx)) '\n']);
k = k+1;
S_k = [S_k A_idx(best_idx)];
A_idx(best_idx)=[];