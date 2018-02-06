function [U] = stage1()
    A = [   
            1 1 1; 
            2 3 1;
            1 -1 -2;
        ];
    
    [m,n] = size(A);
    
    % Verify whether A is a square matrix.
    if m ~= n
        disp("Error: non-square matrix");
        return;
    end
    
     for i = 1:n-1
           for j = i+1:n
                % The pivot is the element in the diagonal line.
                pivot = A(i,i);
                
                % If the pivot is not 1 then simplify the row where the
                % pivot is 1. Any changes made is done to the pivot is done
                % to the other elements in the same row as the pivot.
                if A(i,i) ~= 1
                    A(i,:) = A(i,:)/pivot;
                    pivot = A(i,i);
                end
                
                % Row j = Row j - multiple of row i which is a factor of
                % Row j.
                A(j,:) = A(j,:) - A(i,:) * A(j,i)/pivot;
           end
    end
    
    idx=A==0;
    rank=sum(idx(:));
    
    disp("Rank is of A is: " + rank);
    
    % From the definition of rank in week 8 lecture notes: if rank = m then
    % full rank; if rank < m then rank deficient.
    if rank == m
        disp("Matrix A has full rank.");
    elseif rank < m
        disp("Matrix A is rank deficient.");
    end
    
    U = A;
end