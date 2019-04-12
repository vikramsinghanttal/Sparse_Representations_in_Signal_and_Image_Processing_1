clear all;
clc;
%% setting the parameters
rng(10);
m=150; n=50;                    % dimensions of Grassmanian Matrices
Iter= 1e4;                      % No. of Iterations
dd1	= 0.95;
dd2	= 0.95;

%% Initializing the A, minimum \mu
A   = randn(n,m);
An  = A*diag(1./sqrt(diag(A'*A)));
G   = An'*An;
mu  = sqrt((m-n)/(n*(m-1)));
% An= A./vecnorm(A); 
% G = An'*An;
%% Auxilliary Variables
Res = zeros(Iter,4);

%% Algorithm Implementation
for it=1:Iter
    abs_G_vec=abs(G(:));
    qq = sort(abs_G_vec);
    pos= find(abs_G_vec>qq(round(m*(m-1)*dd1)) & abs(G(:)-1)> 1e-6);
    G(pos)=dd2*G(pos);
    [U,S,V] = svd(G);
    S(n+1:m,n+1:m) = 0;
    G = U*S*V';
    norm_diag = diag(1./sqrt(diag(G)));
    G = norm_diag*G*norm_diag;  
    
    abs_G_pos = abs(G(pos));
    Res(it,:)=[it,mu,mean(abs_G_pos),max(abs_G_pos)];
    fprintf(1,'%6i %12.8f %12.8f %12.8f\n',...
        [Res(it,:)]);  
end
[U,S,V] = svd(G);
Aout = sqrtm(S(1:n,1:n))*U(:,1:n)';
x = 1:1:Iter;
%% Results 
p = figure(1); clf;
set(p,'Position',[415 100 400 200]);
set(gca,'xscale','log')
hold on;
p = semilogx(x,Res(:,2),'b');
p.LineWidth = 2;
p = semilogx(x,Res(:,3),'r');
p.LineWidth = 2;
p = semilogx(x,Res(:,4),'g');
p.LineWidth = 2;
p=legend('optimum \mu','mean \mu','maximum \mu');
lgd1.FontSize = 14;

q = figure(2); clf;
set(q,'Position',[830 100 400 200]);
hold on;
Gi = (abs(An'*An));
qq = sort(Gi(:));
mu_init = sort(qq(1:end-m));

Gf = (abs(Aout'*Aout));
ss = sort(Gf(:));
mu_fin = sort(ss(1:end-m));

q = plot(mu_init,'b');
q.LineWidth = 2;
q = plot(mu_fin,'r');
q.LineWidth = 2;
q = plot(mu*ones(length(mu_init),1),'g');
q.LineWidth = 2;
grid on;
% axis([0 m*(m-1) 0 0.6]);
legend('Initial \mu','FInal \mu','Optimum \mu');