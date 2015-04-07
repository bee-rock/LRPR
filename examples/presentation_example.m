%% Small example for sanity checking
clear;
clc;
% Small example from the beginning of presentation
% Alpha,beta,gamma,...

i = [ 2 6 3 4 4 5 6 1 1];
j = [ 1 1 2 2 3 3 3 4 6];

n = 6;
G = sparse(i,j,1,n,n);

c = sum(G,1);
k = find(c~=0);

D = sparse(k,k,1./c(k),n,n);

% Transition matrix
Pbar = G*D;

[i,j,s] = find(c==0);
d = sparse(i,j,s,1,n);

% Dangling node vector
u = ones(n,1)/n;

% Teleportation vector
v = ones(n,1)/n;

alpha = 0.85;
I = speye(n,n);
% Solutions: Matlab solution and power method
x0 = (I-alpha*(Pbar + u*d))\((1-alpha)*v);

% Power method with function handle
Af = @(x) alpha*Pbar*x + alpha*u*(d*x) + (1-alpha)*v;
[x1,iters] = LRPR_PowerMethod(Af,ones(n,1)/n,.0001,100);



