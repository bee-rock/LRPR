% The Power Method

function [xk,iters] = LRPR_PowerMethod(A,x0,tol,max_iter)


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
    
    if(norm((yk/muk) - xk) < tol)
        display('Reached convergence')
        display(k)
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

