function [distance, C, D] = svdgup(A, B, r, V_init, options)
%svdgup Finds upper bound of the distance to the nearest skew-symmetric
% matrix pencil C-xD of rank <= r to the skew-symmetric matrix pencil A-xB
% using the SVD and GUPTRI.
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
        A double {mustBeNumeric, mustBeSkewSym(A)};
        B double {mustBeNumeric, mustBeSkewSym(B), mustBeEqualSize(A,B)};
        r double {mustBeNonnegative, mustBeWithin(r,A)};
        V_init double;
        options.max_iter double {mustBeNonnegative} = 3000;
        options.epsilon double {mustBeNonnegative} = 10^-7;
        options.phase_1 logical = false;
        options.epsilon_1 double = 10^-3;
    end

    % Constants
    [n,~] = size(A);
    s = floor(r/2);
    N = options.max_iter;
    term_crit = options.epsilon * norm([A;B],"fro");

    % Permutation matrix (transpose)
    Perm = zeros(s*s,s*s);
    for i=1:s % ith row of blocks
        for j=1:s % jth column of blocks
            % This is the block p(i,j)
            % It should have a 1 in position (j,i) (inside the block)
            Perm((i-1)*s+j, (j-1)*s+i) = 1;
        end
    end
    I = eye(s);

    % Initial V
    V = V_init;

    distance = NaN;
    prevdist = NaN;
    mindist = NaN;
    C_min = zeros(n);
    D_min = zeros(n);
    for i=1:N

        % Solve for W using SVD
        [R,S,T] = svd(V);
        S1 = S(1:s,:);

        % Solve for W0...
        E = R'*A*conj(R);
        E11 = E(1:s,1:s);
        E12 = E(1:s,s+1:n);
        
        % Y = T'*W0.'*conj(R)
        Y1 = zeros(s);
        for k = 1:s
            for l = 1:s
                Y1(k,l) = E11(k,l)/(S1(k,k)+S1(l,l));
            end
        end
        
        Y2 = S1\E12;
        Y = [Y1 Y2];
        W0 = R*Y.'*T.';

        % Solve for W1...
        E = R'*B*conj(R);
        E11 = E(1:s,1:s);
        E12 = E(1:s,s+1:n);
        
        % Y = T'*W1.'*conj(R) (s by n)
        Y1 = zeros(s);
        for k = 1:s
            for l = 1:s
                Y1(k,l) = E11(k,l)/(S1(k,k)+S1(l,l));
            end
        end
        
        Y2 = S1\E12;
        Y = [Y1 Y2];
        W1 = R*Y.'*T.';
    
        % Solve for V using GUPTRI
        [S,T,P,Q,~] = pguptri(W0,W1,zeros=true); 

        S1 = S(1:s,:);
        T1 = T(1:s,:);
        T2 = T(s+1:2*s,:);
    
        S23 = S(s+1:n,:);
        T3 = T(2*s+1:n,:);

        % S2, S3 and T3 should be zero. If not, something is wrong.
        assert(isequal(S23, zeros(n-s,s)), "S2-3 non-zero")
        assert(isequal(T3, zeros(n-2*s,s)), "T3 non-zero")
    
        % Calculate E and F
        E = P'*A*conj(P);
        F = P'*B*conj(P);
        
        % Blocks of E and F
        E11 = E(1:s, 1:s);
        E13 = E(1:s, 2*s+1:n);
        E21 = E(s+1:2*s, 1:s);
        F11 = F(1:s, 1:s);
        F12 = F(1:s, s+1:2*s);
        F22 = F(s+1:2*s, s+1:2*s);
        F13 = F(1:s, 2*s+1:n);
        F23 = F(s+1:2*s, 2*s+1:n);
    
        % Solve for X = [X1; X2; X3]; X = P'*V*conj(Q)
        X3 = -([S1; T1; T2]\[E13; F13; F23]).';
        
        % For X1 and X2 we will be using vec...
    
        C1 = kron(S1,I);
        RHS1 = reshape(E21,[],1); % This is vec(E21)
        C2 = (kron(T2,I) - kron(I,T2)*Perm);
        RHS2 = reshape(F22,[],1);
    
        X2vec = [C1; C2]\[RHS1; RHS2];
    
        C3 = (kron(S1,I) - kron(I,S1)*Perm);
        RHS3 = reshape(E11,[],1);
        C4 = (kron(T1,I) - kron(I,T1)*Perm);
        RHS4 = reshape(F11,[],1);
        C5 = kron(T2,I);
        RHS5 = reshape(F12,[],1) + kron(I,T1)*Perm*X2vec;
        
        X1vec = [C3; C4; C5]\[RHS3; RHS4; RHS5];
    
        X1 = reshape(X1vec,s,s);
        X2 = reshape(X2vec,s,s);
    
        X = [X1; X2; X3];
        V = P*X*Q.';

        % Calc distance
        distance = norm([A;B] - [V*W0.'-W0*V.'; V*W1.'-W1*V.'], "fro");

        % In first iteration and whenever distance is less than recorded
        % minimum, update minimum
        if i == 1 || distance < mindist
            mindist = distance;
            C_min = V*W0.'-W0*V.';
            D_min = V*W1.'-W1*V.';
        end

        if isnan(distance)
            break;
        end

        % Check termination criterion
        if abs(prevdist - distance) < term_crit
            break;
        end

        prevdist = distance;
    end

    % return best solution found
    C = C_min;
    D = D_min;
    distance = mindist;

end

function mustBeSkewSym(M)
    if ~isequal(M,-M.')
        error("Matrices must be skew symmetric.")
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