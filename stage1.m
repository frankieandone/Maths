% U1462480 Farzad Merzadyan
function [U] = stage1(A)
    
    % A = [ 1 1 1; 2 3 1; 1 -1 -2];
    
    [m,n] = size(A);
    
    % Verify whether A is a square matrix.
    if m ~= n
        error("Error: non-square matrix.");
        return;
    end
    
    % A is a square matrix then m=n therefore m and n can be used
    % interchangeably.
    for i = 1:n-1
        % The pivot is the element in the diagonal line.
        pivot = A(i,i);
        
        % If the pivot is not 1 then simplify the row where the
        % pivot is 1. Any changes made is done to the pivot is done
        % to the other elements in the same row as the pivot.
        if pivot ~= 1
            A(i,:) = A(i,:)/pivot;
            pivot = A(i,i);
        end
        
        for j = i+1:n
            % Row j = Row j - multiple of row i which is a factor of
            % Row j.
            A(j,:) = A(j,:) - A(i,:) * A(j,i)/pivot;
        end
    end
    
    % U represents upper echelon form of matrix A.
    U = A;
end