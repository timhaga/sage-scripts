#!/usr/bin/env sage

def v_to_matrix(v, N):
    return matrix([list(v)[i:i+N+1] for i in xrange(0, len(list(v)), N+1)])

def Eigensteps(M):
    N = M.ncols()
    return matrix([
        sorted((M[:,0:k]*M[:,0:k].transpose()).eigenvalues())
        for k in range(N+1)]).transpose()

def Psi(N, d, v):
    w = matrix(N-d, N+1)
    for i in xrange(1,N-d+1):
        for n in xrange(d+i,N+1):
            w[N-d-i, n] = N
        for n in xrange(i,d+i):
            w[N-d-i, n] = v[n-i,N-n]
    return w

load("PolytopeOfEigensteps.sage")
load("FUNTFFromEigensteps.sage")

N = 7
d = 2

# Create polytope of FUNTF eigensteps
P = PolytopeOfEigensteps(N,d,d)
print "Polytope of FUNTF eigensteps:"
print P

for v in P.vertices():
    v = v_to_matrix(v, N)
    #print "Some point in P:"
    #print v
    #print "Psi_{%d,%d}(v):" % (N, d)
    #print Psi(N, d, v)

    # Create FUNTF
    F = FUNTFFromEigensteps(v/d)
    #print "An associated FUNTF:"
    #print F
    assert(d*Eigensteps(F) == v)

    # Corresponding parseval frame
    F_p = F * sqrt(d/N)

    # Extend to ONB using Gram Schmidt
    M = F_p
    for i in xrange(0, N):
        U = span(M)
        e_i = matrix(1,N)
        e_i[0,i] = 1
        e_i = e_i[0]
        
        if e_i not in U:
            x = e_i
            for j in xrange(0, M.nrows()):
                x -= (e_i * M[j]) * M[j]
            x = x / x.norm()
            M = M.stack(x)

    assert(M.is_unitary())

    # Obtain Naimark complement
    G_p = M[d:N,:]

    # Correspnding FUNTF
    G = G_p * sqrt(N/(N-d))
    #print "A Naimark complementary FUNTF:"
    #print G
    #print "Eigensteps of reversed complement:"
    #print (N-d)*Eigensteps(G[:,::-1])
    print "Conj. holds:", Psi(N, d, v) == (N-d)*Eigensteps(G[:,::-1])
