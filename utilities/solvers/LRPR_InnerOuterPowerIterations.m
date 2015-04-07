% MATLAB prototype for the inner-outer algorithm from
% "An inner-outer iteration for computing pagerank"
%
% See innout-web for an efficient implementation by authors


function xk = LRPR_InnerOuterPowerIterations(P,alpha,beta,tau,eta,v,Pbar,d,u,max_iter)
    Af = @(x) alpha*Pbar*x + alpha*u*(d*x) + (1-alpha)*v;
    xk = v;
    if isa(P,'function_handle')
        yk = P(xk);
    else
        yk = P*xk;
    end

while norm(alpha*yk + (1-alpha)*v - xk,1) >= tau
    f = (alpha - beta)*yk + (1-alpha)*v;
    % Require at least 1 step in the inner iteration
    xk = f + beta*yk;
    
    if isa(P,'function_handle')
        yk = P(xk);
    else
        yk = P*xk;
    end
    counter = 1;
    while norm(f + beta*yk - xk,1) >= eta
        counter = counter + 1;
        xk = f + beta*yk;
        if isa(P,'function_handle')
            yk = P(xk);
        else
            yk = P*xk;
        end
    end
    % ie: Im = 1 from Algorithm
    if counter == 1
        [xk,~] = LRPR_PowerMethod_pagerank(Af,alpha*yk + (1-alpha)*v,tau,max_iter);
    end
end

xk = alpha*yk + (1-alpha)*v;

