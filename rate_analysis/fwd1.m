function [log_prob] = fwd1(y, m, f, A, sigma, s_0)

Nsamples = length(y);
Nstates = length(f);
F = -Inf*ones(Nstates,Nsamples); 
F_prev = -Inf*ones(Nstates,Nsamples); 

for n = 1:Nsamples  
    for s = find(A(s_0,:)==1)
        F(s,n) = -sum((y(1:n)-f(s)).^2)/(2*sigma^2);
    end
end
F_prev = F;

for ell = 2:m   
    for n = ell:Nsamples
        for s = 1:Nstates 
            F(s,n) = F(s,n-1);
            for s_prev = find(A(:,s)==1)'
                F(s,n) = elnsum(F(s,n), F_prev(s_prev,n-1)); 
            end
            F(s,n) = F(s,n) - (y(n)-f(s))^2/(2*sigma^2);
        end
    end
    F_prev = F;
    ell
end

log_prob = -Inf;
for s = 1:Nstates
    log_prob = elnsum(log_prob,F(s,length(y)));
end

end