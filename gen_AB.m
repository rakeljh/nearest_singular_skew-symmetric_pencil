function [A,B] = gen_AB(n, options)
%gen_AB generates two full rank skew-symmetric complex matrices of size n
    arguments
        n int64
        options.complex logical = true
    end
    
    assert(mod(n,2)==0, "n must be even")
    
    A = zeros(n);
    B = zeros(n);
    
    % Redo if not full rank
    while rank(A) ~= n && rank(B) ~= n
        % Generate real random matrices sampling in [-1,1]
        A = 2*rand(n)-1;
        B = 2*rand(n)-1;

        % Add complex part (sampling in [-i,i]) if requested
        if options.complex
            A = A + 1i*(2*rand(n)-1);
            B = B + 1i*(2*rand(n)-1);
        end
        
        % Make A and B skew symmetric
        A = (A - A.')/2;
        B = (B - B.')/2;
    end

end

