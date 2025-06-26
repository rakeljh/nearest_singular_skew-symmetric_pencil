function [distance, C, D] = vecvec(A, B, r, V_init, options)
%vecvec Finds nearest skew-symmetric matrix pencil C-xD of rank <= r to
% the matrix pencil A-xB, using the "vec-trick".
%
%   Input:
%       A                   :   First coefficient of matrix pencil A-xB.
%       B                   :   Second coefficient of matrix pencil A-xB.
%       r                   :   Rank of matrix to search for, in [2,n).
%       V_init              :   Initial guess for the n by s matrix V.
%       max_iter (optional) :   Maximum number of iterations, default 3000.
%       epsilon (optional)  :   Sensitivity of termination criteria,
%                               default 10^-7.
%
%   Output:
%       distance    :   Distance between matrix pencils A-xB and C-xD, in
%                       the sense of (|A-C|^2 + |B-D|^2)^(1/2) using the
%                       Frobenius norm.
%       C           :   First coefficient of matrix pencil C-xD.
%       D           :   Second coefficient of matrix pencil C-xD.

    arguments
        A double {mustBeNumeric, mustBeSquare(A)};
        B double {mustBeNumeric, mustBeSquare(B), mustBeEqualSize(A,B)};
        r double {mustBeNonnegative, mustBeWithin(r,A)};
        V_init double;
        options.max_iter double {mustBeNonnegative} = 3000
        options.epsilon double {mustBeNonnegative} = 10^-7;
        options.phase_1 logical = false;
        options.epsilon_1 double = 10^-3;
    end

    % Constants
    [n,~] = size(A);
    s = floor(r/2);
    N = options.max_iter;
    term_crit = options.epsilon * norm([A;B],"fro");
    if options.phase_1
        phase = 1;
    else
        phase = 2;
    end
    term_crit_ph1 = options.epsilon_1 * norm([A;B],"fro");

    % Initial V
    V = rand(n,s,like=[A;B]); 
    % V should be real if both A and B are, otherwise complex.

    % Vectorize
    A_vec = reshape(A,[],1);
    B_vec = reshape(B,[],1);

    % Permutation matrix (transpose)
    P = zeros(n*s,n*s);
    for i=1:n % ith row of blocks
        for j=1:s % jth column of blocks
            % This is the block p(i,j)
            % It should have a 1 in position (j,i) (inside the block)
            % The block contains P((i-1)*s+1:i*s, (j-1)*n+1:j*n)
            % The jth row of the block is the row (i-1)*s+j of P
            % and the ith column is (j-1)*n+i
            P((i-1)*s+j, (j-1)*n+i) = 1;
        end
    end
    
    I = eye(n);
    distance = NaN;
    prevdist = NaN;
    for i=1:N

        if phase == 1

            % Approximate V, solve for W1, approximate V, solve for W0, ...
            W0 = svd_solver(V,A);
            V = svd_solver(W0,-A);
            W1 = svd_solver(V,B);    
            V = svd_solver(W1,-B);

        else
            % Solve for W
            V_tilde = kron(I, V)*P - kron(V, I);
        
            dV = decomposition(V_tilde); % saves time
            W0_vec = dV\A_vec;
            W1_vec = dV\B_vec;
        
            W0 = reshape(W0_vec,n,s);
            W1 = reshape(W1_vec,n,s);
        
            % Solve for V
            W0_tilde = kron(W0, I) - kron(I, W0)*P;
            W1_tilde = kron(W1, I) - kron(I, W1)*P;
        
            V_vec = [W0_tilde; W1_tilde]\[A_vec; B_vec];
        
            V = reshape(V_vec,n,s);
        end
    
        % Calc distance
        distance = norm([A;B] - [V*W0.'-W0*V.'; V*W1.'-W1*V.'], "fro");

        % Check termination criterion
        if abs(prevdist - distance) < term_crit
            break;
        elseif phase == 1 && abs(prevdist - distance) < term_crit_ph1
            phase = 2;
        end

        prevdist = distance;
    end


    C = V*W0.'-W0*V.';
    D = V*W1.'-W1*V.';
end

function mustBeSquare(M)
    [m,n] = size(M);
    if m ~= n
        error("Matrices must be square.")
    end
end

function mustBeEqualSize(M,N)
    if ~isequal(size(M),size(N))
        error("Matrices must be of equal size.")
    end
    if any([size(M),size(N)] <= 2)
        error("Size of input pencil must be greater than 2 by 2.")
    end
end

function mustBeWithin(r,M)
    if r < 2 || r >= length(M)
        error("Rank r must be in [2,n).")
    end
end