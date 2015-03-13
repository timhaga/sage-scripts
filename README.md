# a collection of sage scripts

## frames.sage — calculations with frames and eigensteps

`frames.sage` can do the following things:

 - `Eigensteps(F)` calculates the eigensteps of a frame F
 - `EigenstepsGram(F)` calculate outer eigensteps
 - `FrameFromEigensteps(lam)` constructs a frame from a given
   sequence of eigensteps
 - `PolytopeOfEigensteps(mu, lam)` constructs the polytope of
   eigensteps for given sequences of norm squares and eigenvalues
 - `PolytopeOfEigenstepsENT(N, d, mu)` does the same for equal
   norm tight frames
 - `anInteriorPointENT(N, d, mu)` calculates an interior point
   of a polytope of eigensteps of equal norm tight frames
 - `Psi(lam)` and `Phi(lam)` apply affine isomorphisms mentioned
   in our paper "Polytopes of Eigensteps of Finite Equal Norm
   Tight Frames"
 - `v_to_matrix(v, N)` and `matrix_to_v(lam)` convert between
   different formats of sequences of eigensteps (matrices can't
   be points of polytopes in sage)
 - `NaimarkComplement(F)` calculates a Naimark complement of
   a Bessel sequence
 - `CanonicalCompletion(F)` calculates the canonical completion
   of a Bessel sequence to a tight frame
