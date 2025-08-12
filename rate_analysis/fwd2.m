function [log_prob] = fwd2(y, m, f, sigma, s_vec)

Nsamples = length(y);
F = -Inf*ones(1,Nsamples); 
F_prev = -Inf*ones(1,Nsamples); 

for n = 1:Nsamples  
    F(n) = -sum((y(1:n)-f(s_vec(1))).^2)/(2*sigma^2);
end
F_prev = F;

for ell = 2:m   
    for n = ell:Nsamples
        F(n) = elnsum(F(n-1), F_prev(n-1)) - (y(n)-f(s_vec(ell)))^2/(2*sigma^2); 
    end
    F_prev = F;
end

log_prob = F(length(y));

end