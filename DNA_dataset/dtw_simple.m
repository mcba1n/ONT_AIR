function [best_cost, dtw_matrix] = dtw_simple(x,y)
% simple dtw algorithm from wikipedia

m=length(x);
n=length(y);

dtw_matrix = Inf*ones(n+1,m+1);
dtw_matrix(1,1) = 0;

for j = 2:m+1
   for i = 2:n+1
      dtw_matrix(i,j) =  (y(i-1) - x(j-1))^2 + min([...
                          dtw_matrix(i-1,j), ... %insertion
                          %dtw_matrix(i,j-1), ... %deletion
                          dtw_matrix(i-1,j-1)]); %match
   end
end

best_cost = dtw_matrix(n+1,m+1);

end

