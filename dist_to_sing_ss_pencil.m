function [distance, C, D] = dist_to_sing_ss_pencil(A, B, r, options)
%dist_to_sing_ss_pencil Finds nearest skew-symmetric matrix pencil C-xD of
% rank <= r, with r in [2,n), to the n by n matrix pencil A-xB, using the
% Frobenius norm.
%
%   Input:
%       A                   :   First coefficient of matrix pencil A-xB.
%       B                   :   Second coefficient of matrix pencil A-xB.
%       r                   :   Rank of matrix to search for, in [2,n).
%       max_iter (optional) :   Maximum number of iterations, default 3000.
%       epsilon (optional)  :   Sensitivity of termination criteria,
%                               default 10^-7.
%       u_bound (optional)  :   If true, returns an upper bound in less
%                               time. Requires GUPTRI in MATLAB path.
%                               Default false. Recommended only for pencils
%                               of size 20 by 20 or larger.
%       phase_1 (optional)  :   If true, improves the starting guess by
%                               cheap iterations before entering the main
%                               algorithm. ONLY for skew-symmetric pencils.
%                               Default false.
%       epsilon_1 (optional):   Sensitivity of termination criteria for
%                               phase 1. Default 10^-3.
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
        r double {mustBeNonnegative, mustBeWithin(r,A)} = 0;
        options.max_iter double {mustBeNonnegative} = 3000;
        options.epsilon double {mustBeNonnegative} = 10^-7;
        options.u_bound logical = false;
        options.phase_1 logical = false;
        options.epsilon_1 double = 10^-3;
    end

    % Unpack options
    N = options.max_iter;
    eps = options.epsilon;
    GUPTRI = options.u_bound;
    ph_1 = options.phase_1;

    n = length(A);
    s = floor(r/2);
    V = rand(n,s,like=[A;B]);

    if ph_1
        if isequal(A,-A.') && isequal(B,-B.')
            term_crit_ph1 = options.epsilon_1 * norm([A;B],"fro");

            prevdist = NaN;
            for i=1:100 % max iter for phase 1 is 100
                W0 = svd_solver(V,A);
                V = svd_solver(W0,-A);
                W1 = svd_solver(V,B);    
                V = svd_solver(W1,-B);

                dist = norm([A;B] - [V*W0.'-W0*V.'; V*W1.'-W1*V.'], "fro");

                if abs(prevdist - dist) < term_crit_ph1
                    break;
                end

                prevdist = dist;
            end
        else
            warning("Phase 1 not compatible with non skew-symmetric" + ...
                "pencils, phase 1 not used.");
        end
    end

    % Check skew symmetry and choose method.
    if isequal(A,-A.') && isequal(B,-B.')
        if GUPTRI
            warning("off") % Produces warnings about rank deficiency.
            [distance, C, D] = svdgup(A, B, r, V, ...
                max_iter=N, epsilon=eps);
            warning("on")
        else
            [distance, C, D] = svdvec(A, B, r, V, ...
                max_iter=N, epsilon=eps);
        end
    else
        warning("Input matrices not skew symmetric, " + ...
            "performance will be affected.")
        warning("off") % Produces warnings about rank deficiency.
        [distance, C, D] = vecvec(A, B, r, V,...
            max_iter=N, epsilon=eps);
        warning("on")
    end
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
    if 2 > r || r >= length(M)
        error("Rank r must be in [2,n).")
    end
end

function X = svd_solver(V,C)
%svd_solver returns X as the solution to VX.' - XV.' = C.
    
    % V and X have the same dim, n by s
    [n,s] = size(V);

    % SVD of V
    [R,S,T] = svd(V);
    S1 = S(1:s,:);
    E = R'*C*conj(R);
    E11 = E(1:s,1:s);
    E12 = E(1:s,s+1:n);
    
    % Solve for Y = T'*X.'*conj(R) (s by n)
    Y1 = zeros(s);
    for k = 1:s
        for l = 1:s
            Y1(k,l) = E11(k,l)/(S1(k,k)+S1(l,l));
        end
    end
    
    Y2 = S1\E12;
    
    Y = [Y1 Y2];
    
    % Solve for X
    X = R*Y.'*T.';

end