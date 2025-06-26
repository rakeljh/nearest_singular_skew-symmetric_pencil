% Choose size and rank
n = 10;
r = n-1;

% Generate input pencil
[A,B] = gen_AB(n);

% Solve
[dist, C, D] = dist_to_sing_ss_pencil(A,B,r);