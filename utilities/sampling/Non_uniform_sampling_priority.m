% Non uniform Sampling Scheme as outlined
% in Fast PageRank approximation by adaptive sampling
% by Wenting Liu · Guangxia Li · James Cheng
% 
%
% Input: Square sparse matrix A, Initialized priority queue
% 
% Author: Brock Hargreaves
%
% This version uses an efficient priority queue from a mex compiled 
% cpp implementation.

function Atilde = Non_uniform_sampling_priority(A,pq)
    n = size(A,1);
    N = nnz(A);
    s = floor(N/2);
    theta = (8*log(n))^2/sqrt(n);
    Z = 0;
    
    [i,j,A_ij] = find(A);
   
    for l = 1:N
        Z = Z + A_ij(l)^2;
        r_ij = rand(1,1);
        mu = s*A_ij(l)^2;
        k_ij = max(mu/r_ij,(mu/r_ij^2)*theta^2);
        
        % Since we are using a Max heap implementation of the queue, we
        % push on the negative key and remove all those entries which are
        % greater than -Z, rather than remove all those entries smaller
        % than Z
        pq_push(pq, l, -k_ij);
        
        [~,cost] = pq_top(pq);
        while cost > -Z
            pq_pop(pq);
            [~,cost] = pq_top(pq);
        end
        
    end
    
    newsize = pq_size(pq);
    idx = zeros(newsize,1);
    
    for k = 1:newsize
        [idx(k),~] = pq_pop(pq);
    end
    
    i_new = i(idx);
    j_new = j(idx);
    Atilde = sparse(i_new,j_new,A_ij(idx),n,n);

end  
    
    
    