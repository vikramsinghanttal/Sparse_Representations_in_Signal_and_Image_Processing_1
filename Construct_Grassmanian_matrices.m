clear all;
clc;
%% setting the parameters
m=150; n=50;                    % dimensions of Grassmanian Matrices
Iter= 1e4;                      % No. of Iterations
dd1	= 0.95;
dd2	= 0.95;

%% Initializing the A, minimum \mu
A   = randn(n,m);
G   = A'*A;
G   = G./diag(G);
mu  = sqrt((m-n)/(n*(m-1)));

A = [16, -2, 15, 3;
     5,   6,  8, 8;
     9,   4, 11,11;
     4,   12,10, 1];
 An= A./vecnorm(A); 
 G = An'*An;