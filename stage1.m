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
    % n represents the outer boundaries (since m=n) of the matrix.
    for i = 1:n-1
        % A(i:n,i) gets all elements from row i to n in column i.
        % abs(A(i:n,i)) gets the absolute values for these elements.
        % max(abs(A(i:n,i))) gets the largest absolute value.
        % v (~) holds the largest absolute value. It is not used therefore
        % ~ is used as an anonymous variable; v is needed to make sure that
        % return of max(abs(A(i:n,i))) is dimensionally compatible.
        % k represents position of v in relation to row i. e.g.: previous
        % row to i, next row to i, etc.
        [~, k] = max(abs(A(i:n,i)));
        
        % If k = 1 then it means same row, no point checking the row in
        % that case. k > 1 prevents useless iteration.
        if k > 1
            % Variable temp holds a temporary copy of current row.
            temp = A(i,:);
            % Move next row into current row.
            A(i,:) = A(i+k-1,:);
            % Copy the temporary row into the next row completing the swap.
            A(i+k-1,:) = temp;
        end
        
        % The pivot is the element in the diagonal line.
        pivot = A(i,i);
        
        % Guard against divide by zero errors.
        if pivot == 0
            break;
        end
        
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
    
    % tnZ represents total non-zero rows in U.
    % Assume matrix has full rank.
    tNZ = n;
    for i = 1:n
        % isZero is true if all elements in row i are zero values.
        isZeroRow = all(A(i,:) == 0);
        % If row i are all zero values then decrement rank by 1.
        if isZeroRow == 1 && tNZ > 0
            tNZ = tNZ - 1;
        end
    end
    
    % Using definitions from lecture notes from week 8: 
    % full rank is when rank = n.
    % rank deficient is when rank < n.
    if tNZ == n
        disp("Full rank. Rank: " + tNZ);
    elseif tNZ < n
        disp("Rank deficient. Rank: " + tNZ);
    elseif tNZ < 0
        error("Error: rank < 0. Rank: " + tNZ);
    end
end