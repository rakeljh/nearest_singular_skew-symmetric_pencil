# nearest_singular_skew-symmetric_pencil
Find nearest skew-symmetric matrix pencil of rank &lt;=r. For an n by n input pencil, r must be in [2,n).

Example usage found in example.m. The file gen_AB.m is used to generate a random input pencil for testing. The main function file is dist_to_sing_ss_pencil, which uses vecvec, svdvec or svdgup depending on input.

Note that svdgup (u_bound option in dist_to_sing_ss_pencil) requires the MCS Toolbox https://www.umu.se/forskning/projekt/stratigraph-and-mcs-toolbox/.
