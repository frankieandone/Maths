% U1462480 Farzad Merzadyan
function [U] = stage1(A)
    
    % A = [ 1 1 1; 2 3 1; 1 -1 -2];
    
    [m,n] = size(A);
    
    % Verify whether A is a square matrix.
    if m ~= n
        error("Error: non-square matrix.");
        return;
    end
    
    for i = 1:n-1
           for j = i+1:n
                % The pivot is the element in the diagonal line.
                pivot = A(i,i);
                
                % If the pivot is not 1 then simplify the row where the
                % pivot is 1. Any changes made is done to the pivot is done
                % to the other elements in the same row as the pivot.
                if pivot ~= 1
                    A(i,:) = A(i,:)/pivot;
                    pivot = A(i,i);
                end
                
                % Row j = Row j - multiple of row i which is a factor of
                % Row j.
                A(j,:) = A(j,:) - A(i,:) * A(j,i)/pivot;
           end
    end
    
    U = A;
    
    % Count the instances of non-zero rows in matrix A.
    idx=U==0;
    RankOfU=sum(idx(:));
    
    disp("Rank is of matrix U: " + RankOfU);
    
    % From the definition of rank in week 8 lecture notes: if rank = m then
    % full rank; if rank < m then rank deficient.
    if RankOfU == m
        disp("Matrix U has full rank.");
    elseif RankOfU < m
        disp("Matrix U is rank deficient.");
    elseif RankOfU > m
        % Catches improbable events when RankOfU > m.
        % Should be impossible but catching just in case.
        error("Error: rank cannot be greater than m.");
    end
end