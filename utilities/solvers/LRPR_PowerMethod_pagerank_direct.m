% The Power Method for computing Pagerank with criterion for direct
% sampling

function [xk,iters] = LRPR_PowerMethod_pagerank_direct(A,x0,tau,max_iter)

%Debugging
verbose = 0;

if norm(x0,1) > 1 + 10e-12
    display('Warning: Initial guess is not normalized');
end

xk = x0;

for k = 1:max_iter
    if isa(A,'function_handle')
        yk = A(xk);
    else
        yk = A*xk;
    end

    muk = norm(yk,1);
    
    if(norm(xk - yk/muk,1) < tau)
        if verbose
            display('Reached convergence')
            display(k)
        end
        iters = k;
        xk = yk/muk;
        return
    else if k == max_iter
            iters = max_iter;
            display('Maximum iterations')
        end
    end
    xk = yk/muk;
end

