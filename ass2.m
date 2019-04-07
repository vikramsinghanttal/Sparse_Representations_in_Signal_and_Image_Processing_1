clear all;
clc
m=4;
A = [16, -2, 15, 13;
     5,   6,  8, 8;
     9,   4, 11,12;
     4,   12,10, 1];
 A  = A./vecnorm(A);
 A  = A./vecnorm(A);
G  = (A'*A);


Gsorted = sort(G,2,'descend');
Gsorted(:,1)=[];
Gcumsum = cumsum(Gsorted,2);
max_Gcumsum = max(Gcumsum);