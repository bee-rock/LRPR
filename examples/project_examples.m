%% Testing all methods
%
% Each example is a full run of this file and then saved in a matlab
% workspace.
% Only two lines need to be changed in this file to test various examples
% 1. Working directory: wd
% 2. example
clear;
clc;

% Set working directory 
wd = '/Users/brockhargreaves/Documents/MATLAB/LRPR';
example = 'ubc-cs';



fprintf('Setting experiment parameters.\n')
% Damping parameter
alpha = 0.85;

% IOPI (Inner outer Power Iteration) Parameters
beta  = 0.5;
eta = 10e-4;

% Convergence criterion
tau = 10e-5;

% Read in data and construct relevant sparse matrices
example_path = strcat('/utilities/innout-small/data/', example);
G = bvgraph(strcat(wd, example_path));
G = sparse(G);

c = sum(G,1);
k = find(c~=0);
n = size(G,1);
D = sparse(k,k,1./c(k),n,n);

% Transition matrix
Pbar = G*D;
 
% Create dangling and teleportation probability vectors
[i,j,s] = find(c==0);
d = sparse(i,j,s,1,n);

% Dangling node vector
u = ones(n,1)/n;

% Teleportation vector
v = ones(n,1)/n;

I = speye(n,n);

% Begin Experiment
% Power method
Af = @(x) alpha*Pbar*x + alpha*u*(d*x) + (1-alpha)*v;
 

tic;
fprintf('Solving example with power method.\n')
[x_PI,iters] = LRPR_PowerMethod_pagerank(Af,ones(n,1)/n,tau,10000);
toc;


% Sample the Transition matrix

% Note: If you decide to run this multiple times, you may need to clear the
% MATLAB workspace, otherwise MATLAB might crash. It's a weird memory bug,
% though pq_delete(pq) should be handling it.
N = nnz(Pbar);
pq = pq_create( N ); 
tic
fprintf('Sampling the transiiton matrix. This may take some time.\n')
PbarTilde = Non_uniform_sampling_priority(Pbar,pq);
toc
pq_delete(pq);


% Power method with direct sampling
Aft = @(x) alpha*PbarTilde*x + alpha*u*(d*x) + (1-alpha)*v;
 
tic;
fprintf('Solving example with Power method and a subsampled transition matrix.\n')
[x_DSPI,iters] = LRPR_PowerMethod_pagerank_direct(Aft,ones(n,1)/n,tau,10000);
toc;

% Inner Outer

P = @(x) Pbar*x + u*(d*x);
tic;
fprintf('Solving example with inner outer iterations.\n')
x_IO = LRPR_InnerOuter(P,alpha,beta,tau,eta,v);
toc;


% Inner Outer with direct sampling

P = @(x) PbarTilde*x + u*(d*x);
tic;
fprintf('Solving example with inner outer iterations with subsampled transition matrix.\n')
x_DSIO = LRPR_InnerOuter(P,alpha,beta,tau,eta,v);
toc;

P = @(x) Pbar*x + u*(d*x);
tic
fprintf('Solving example with inner outer power iterations. (Not included in results, not working so well).\n')
x4 = LRPR_InnerOuterPowerIterations(P,alpha,beta,tau,eta,v,Pbar,d,u,10000);
toc;

fprintf('Compute Relative errors using the Power Method as a reference.\n')

rel_err_IO   = norm(x_IO - x_PI,1)/norm(x_PI,1);
rel_err_DSIO = norm(x_DSIO - x_PI,1)/norm(x_PI,1);
rel_err_DSPI = norm(x_DSPI - x_PI,1)/norm(x_PI,1);

fprintf('Computing Spearman rank coefficients.\n')

[rank_vector_PI]   = tiedrank(x_PI);
[rank_vector_IO]   = tiedrank(x_IO);
[rank_vector_DSIO] = tiedrank(x_DSIO);
[rank_vector_DSPI] = tiedrank(x_DSPI);

rho_IO   = sum((rank_vector_IO - mean(rank_vector_IO)).*(rank_vector_PI - mean(rank_vector_PI)))/sqrt(sum((rank_vector_IO - mean(rank_vector_IO)).^2) * sum((rank_vector_PI - mean(rank_vector_PI)).^2));
rho_DSIO = sum((rank_vector_DSIO - mean(rank_vector_DSIO)).*(rank_vector_PI - mean(rank_vector_PI)))/sqrt(sum((rank_vector_DSIO - mean(rank_vector_DSIO)).^2) * sum((rank_vector_PI - mean(rank_vector_PI)).^2));
rho_DSPI = sum((rank_vector_DSPI - mean(rank_vector_DSPI)).*(rank_vector_PI - mean(rank_vector_PI)))/sqrt(sum((rank_vector_DSPI - mean(rank_vector_DSPI)).^2) * sum((rank_vector_PI - mean(rank_vector_PI)).^2));

fprintf('Finished!\n')