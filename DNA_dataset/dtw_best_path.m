function [x,ix] = dtw_best_path(xn,dtw_matrix)
% find the best path in dtw matrix, 
% and the corresponding de-warped signals

i = size(dtw_matrix,1);
j = size(dtw_matrix,2);

ix = [length(xn)];
x = [xn(end)];

while i>2
    % get operation with lowest cost
    [~,I]=min([dtw_matrix(i-1,j),dtw_matrix(i-1,j-1)]);
    %min([dtw_matrix(i-1,j), dtw_matrix(i,j-1), dtw_matrix(i-1,j-1)]);
    
    if I == 1
        % insertion
        i = i-1;
        j = j;
        x = [xn(j-1), x];
    elseif I==2
        % match
        i = i-1;
        j = j-1;
        x = [xn(j-1), x];
    end

    ix = [j-1,ix];
    
    %     elseif I==2
%         % deletion
%         i = i;
%         j = j-1;



    
end

end

