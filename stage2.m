% U1462480 Farzad Merzadyan
function [x] = stage2(A,b)
    % A = [ 1, 2, -3; 2, -1, 1; 4, -1, 1 ];
    % b = [ 9; 0; 4 ];
    
    [m,n] = size(A);
    [m2, n2] = size(b);
    
    % Verify whether A is a square matrix.
    if m ~= n
        error("Error: non-square matrix.");
        return;
    end
    
    % Verify whether matrix A and vector b are compatible.
    if n ~= m2
        error("Error: incompatible dimensions between A and b.");
        return;
    end
    
    % Form augmented matrix AugMatrix.
    AugMatrix = [A,b];
    
    % A is a square matrix then m=n therefore m and n can be used
    % interchangeably.
    for i = 1:n-1
           for j = i+1:n
                % The pivot is the element in the diagonal line.
                pivot = AugMatrix(i,i);
                
                % If the pivot is not 1 then simplify the row where the
                % pivot is 1. Any changes made is done to the pivot is done
                % to the other elements in the same row as the pivot.
                if pivot ~= 1
                    AugMatrix(i,:) = AugMatrix(i,:)/pivot;
                    pivot = AugMatrix(i,i);
                end
                
                % Row j = Row j - multiple of row i which is a factor of
                % Row j.
                AugMatrix(j,:) = AugMatrix(j,:) - AugMatrix(i,:) * AugMatrix(j,i)/pivot;
           end
    end
    
    % U represents upper echelon form of the augmented matrix.
    U = AugMatrix;
    
    % Backwards substitution.
    
    % A is a square matrix then m=n therefore m and n can be used
    % interchangeably.
    % Reminder: [m,n] = size(A);
    % Reminder: [m2, n2] = size(b);
    % Create a vector with the same dimensions as vector b.
    % TODO: comment code.
    x(n, 1:n2) = U(n, n+1:n+n2)/U(n,n);
    for i = n-1:-1:1
        x(i, 1:n2) = (U(i, n+1:n+n2) - U(i, i+1:n)*x(i+1:n, 1:n2)) / U(i,i);
    end
end