% MATLAB prototype for the inner-outer algorithm from
% "An inner-outer iteration for computing pagerank"
%
% See innout-web for an efficient implementation by authors


function xk = LRPR_InnerOuter(P,alpha,beta,tau,eta,v)

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
    while norm(f + beta*yk - xk,1) >= eta
        xk = f + beta*yk;
        if isa(P,'function_handle')
            yk = P(xk);
        else
            yk = P*xk;
        end
    end
end

xk = alpha*yk + (1-alpha)*v;

