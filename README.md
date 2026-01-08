# nearest_singular_skew-symmetric_pencil
Given the $n \times n$ matrix pencil $A - \lambda B$, find nearest skew-symmetric matrix pencil $C - \lambda D$ of rank $\leq r$, where $r \in [2,n)$. Distance is defined as

$\text{dist}(A - \lambda B, C - \lambda D) = \sqrt{ \Vert A - C \Vert_F^2 + \Vert B - D \Vert_F^2}$,

where $\Vert \cdot \Vert_F$ denotes the Frobenius norm.

Example usage found in example.m. The file gen_AB.m is used to generate a random input pencil for testing. The main function file is dist_to_sing_ss_pencil, which uses vecvec, svdvec or svdgup depending on input.

Note that svdgup (u_bound option in dist_to_sing_ss_pencil) requires the MCS Toolbox https://www.umu.se/forskning/projekt/stratigraph-and-mcs-toolbox/.

For matrix polynomials, see https://github.com/rakeljh/dist_sing_skew_pol.
