def Eigensteps(F):
    r"""
    Calculate the (outer) eigensteps of a frame.

    INPUT:

    - ``F`` -- a $d\times N$ matrix with the frame vectors as columns

    OUTPUT:

    a $d \times (N+1)$ matrix of the eigensteps of F
    """
    N = F.ncols()
    return matrix([
        sorted((F[:,0:k]*F[:,0:k].conjugate_transpose()).eigenvalues(), reverse=True)
        for k in range(N+1)]).transpose()

def EigenstepsGram(F):
    r"""
    Calculate the inner eigensteps of a frame.

    INPUT:

    - ``F`` -- a $d\times N$ matrix with the frame vectors as columns

    OUTPUT:

    a list of spectra of partial Gram matrices
    """
    N = F.ncols()
    return [
        sorted((F[:,0:k].conjugate_transpose()*F[:,0:k]).eigenvalues(), reverse=True)
        for k in range(1,N+1)]

def FrameFromEigensteps(lam):
    r"""
    Construct a frame from a sequence of eigensteps.

    INPUT:

    - ``lam`` -- a valid sequence of eigensteps

    OUTPUT:

    a frame with the desired eigensteps
    """
    d = lam.nrows()
    N = lam.ncols() - 1

    F = matrix(SR, d, 0)

    mu_1 = sum(lam.column(1))
    U = matrix.identity(d)

    F = F.augment(sqrt(mu_1)*U.column(0))

    for n in xrange(1,N):
        # step B.1
        V = matrix.identity(d)

        # step B.2
        I = Set(range(d))
        J = Set(range(d))

        for m in reversed(xrange(d)):
            if lam[m,n] in {lam[j,n+1] for j in J}:
                m_J = max({j for j in J if lam[j,n+1]==lam[m,n]})
                I = I - {m}
                J = J - {m_J}

        I = { i+1 for i in I } # switch back to indices >=1
        J = { j+1 for j in J }

        R = len(I)

        pi_I = Permutation(list(I) + list(Set([1..d]).difference(I))).inverse()
        pi_J = Permutation(list(J) + list(Set([1..d]).difference(J))).inverse()

        # step B.3
        v = vector(SR, d)
        w = vector(SR, R)

        for m in I:
            v[pi_I(m)-1] = sqrt(-prod((lam[m-1,n] - lam[j-1,n+1]) for j in J)/ \
                                 prod((lam[m-1,n] - lam[i-1,n]) for i in I-{m}))
        for m in J:
           w[pi_J(m)-1] = sqrt(prod((lam[m-1,n+1] - lam[i-1,n]) for i in I)/ \
                               prod((lam[m-1,n+1] - lam[j-1,n+1]) for j in J-{m}))

        # step B.4
        f = U * V * pi_I.inverse().to_matrix() * v

        # step B.5
        W = matrix(SR, R)
        for i in I:
            for j in J:
                W[pi_I(i)-1, pi_J(j)-1] = \
                        (1/(lam[j-1,n+1] - lam[i-1,n])) * \
                        v[pi_I(i)-1] * w[pi_J(j)-1]

        U = U * V * pi_I.inverse().to_matrix() * \
                block_diagonal_matrix(W, identity_matrix(SR, d-R)) * pi_J.to_matrix()
        F = F.augment(f)

    return F

def PolytopeOfEigensteps(mu, lam):
    r"""
    Construct the polytope of eigensteps

    INPUT:

    - ``mu`` -- a list of norm squares
    - ``lam`` -- a non-increasing spectrum majorizing ``mu``

    OUTPUT:

    The polytope of eigensteps
    """
    N = len(mu)
    d = len(lam)

    # Matrix of variables
    x = matrix(SR, d, N+1)
    for n in xrange(0, N+1):
        for j in xrange(0, d):
            x[j,n] = var("x_%d_%d" % (j+1, n))

    # Collect the interlacing conditions
    ineqs = Set()
    for n in xrange(0,N):
        for j in xrange(0,d):
            ineqs += Set([x[j,n+1] - x[j,n]])
        for j in xrange(0,d-1):
            ineqs += Set([x[j,n] - x[j+1,n+1]])

    # Collect equalities
    eqs = Set()
    # First column is zero, last column is lam_i
    for j in xrange(0,d):
        eqs += Set([x[j,0], x[j,N]-lam[d-j-1]])

    # The other columns sum to mu_n
    for n in xrange(1,N+1):
        S = SR(0)
        for j in xrange(0,d):
            S += x[j,n]
        eqs += Set([S-sum(mu[:n])])

    # Transform the inequalities into a list as needed for Polyhedron()
    ineqs_l = []
    for ineq in ineqs:
        ineq_l = [ineq.polynomial(SR).constant_coefficient()]
        for j in xrange(0, d):
            for n in xrange(0, N+1):
                y = x[j,n]
                ineq_l.append(ineq.coefficient(y,1))
        ineqs_l.append(ineq_l)
    # Same for equalities
    eqs_l = []
    for eq in eqs:
        eq_l = [eq.polynomial(SR).constant_coefficient()]
        for j in xrange(0, d):
            for n in xrange(0, N+1):
                y = x[j,n]
                eq_l.append(eq.coefficient(y,1))
        eqs_l.append(eq_l)

    # Create the polytope
    return Polyhedron(eqns=eqs_l, ieqs=ineqs_l, backend='ppl', base_ring=QQ)

def PolytopeOfEigenstepsENT(N, d, mu=1):
    r"""
    Construct the polytope of eigensteps for equal norm tight frames

    INPUT:

    - ``N`` -- length of the frames
    - ``d`` -- dimension of Hilbert space
    - ``mu`` -- (default: 1) norm square

    OUTPUT:

    The polytope of eigensteps
    """
    return PolytopeOfEigensteps([mu]*N, [N*mu/d]*d)

def anInteriorPointENT(N, d, mu=1):
    r"""
    Construct an interioir point of a polytope of eigensteps of equal norm tight frames

    INPUT:

    - ``N`` -- length of the frames
    - ``d`` -- dimension of Hilbert space
    - ``mu`` -- (default: 1) norm square

    OUTPUT:

    The desired interior point
    """
    p = matrix(d, N+1)
    for i in xrange(1,d+1):
        for n in xrange(i, min(N-d+i,N+1)):
            p[i-1,n] = ((d+n-2*i+1)/d)*mu
        for n in xrange(min(N-d+i,N+1),N+1):
            p[i-1,n] = N/d*mu
    return p.list()

def Psi(lam):
    r"""
    The affine isomorphism $\Psi$ for polytopes of eigensteps of equal norm tight frames

    INPUT:

    - ``lam`` -- a sequence of eigensteps

    OUTPUT:

    The image of ``lam`` under $\Psi$
    """
    d = lam.nrows()
    N = lam.ncols() - 1
    lam_ = matrix(N-d, N+1)
    for i in xrange(1,N-d+1):
        for n in xrange(i,min(d+i,N+1)):
            lam_[i-1, n] = lam[d+i-n-1,N-n]
        for n in xrange(min(d+i,N+1),N+1):
            lam_[i-1, n] = N
    return lam_

def Phi(lam):
    r"""
    The affine isomorphism $\Phi$ for polytopes of eigensteps of equal norm tight frames

    INPUT:

    - ``lam`` -- a sequence of eigensteps

    OUTPUT:

    The image of ``lam`` under $\Phi$
    """
    d = lam.nrows()
    N = lam.ncols() - 1
    X = Set(lam.column(N).list())
    assert(len(X)==1)
    return X[0] * ones_matrix(d, N+1) - lam[::-1, ::-1]

def v_to_matrix(v, N):
    r"""
    Convert a point of a polytope of eigensteps to a sequence of eigensteps matrix

    INPUT:

    - ``v`` -- the point to convert
    - ``N`` -- length of frames

    OUTPUT:

    The corresponding sequence of eigensteps matrix
    """
    return matrix([list(v)[i:i+N+1] for i in xrange(0, len(list(v)), N+1)])
def matrix_to_v(lam):
    r"""
    Convert a sequence of eigensteps matrix to a point $(\lambda_{1,0},\dots,\lambda_{d,N})$

    INPUT:

    - ``lam`` -- a sequence of eigensteps

    OUTPUT:

    The corresponding point
    """
    return lam.list()

def NaimarkComplementPF(F):
    r"""
    Compute a Naimark complement of a Parseval frame

    INPUT:

    - ``F`` -- a Parseval frame

    OUTPUT:

    A Naimark complement of F
    """
    N = F.ncols()
    d = F.nrows()

    # Check that F is a Parseval frame
    assert(F*F.conjugate_transpose() == identity_matrix(d, N))

    # Extend to ONB using Gram Schmidt
    M = F
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
    G = M[d:N,:]

    return G

def CanonicalCompletion(F):
    r"""
    Compute the canonical completion of F to a tight frame

    INPUT:

    - ``F`` -- a Bessel sequence

    OUTPUT:

    The canonical completion
    """
    N = F.ncols()
    d = F.nrows()

    FF = F*F.conjugate_transpose()

    # obtain current spectrum
    ev = sorted(FF.eigenvectors_right(), reverse=True)
    # obtain upper frame bound
    B = ev[0][0]
    # current multiplicity
    K = ev[0][2]

    if K==d:
        return F

    # construct list of remaining eigenvalues and eigenvectors
    lamphi = []
    for e in ev[1:]:
        for phi in e[1]:
            lamphi.append( (e[0], phi) )

    # calculate canonical completion to a B-tight frame
    H = matrix(map(lambda (lam, phi) : sqrt(B-lam)*phi, lamphi)).transpose()

    FH = F.augment(H)

    assert(FH*FH.conjugate_transpose() == B*identity_matrix(d))
    
    return FH

def NaimarkComplement(F):
    r"""
    Compute a Naimark complement of a Bessel sequence

    INPUT:

    - ``F`` -- a Bessel sequence

    OUTPUT:

    A Naimark complement of F
    """
    N = F.ncols()
    d = F.nrows()

    F_t = CanonicalCompletion(F)
    N_t = F_t.ncols()

    FF_t = F_t*F_t.conjugate_transpose()
    B = FF_t[0,0]
    assert(FF_t == B*identity_matrix(d))

    M = F_t/sqrt(B)

    # Extend to ONB using Gram Schmidt
    for i in xrange(0, N_t):
        U = span(M)
        e_i = matrix(1,N_t)
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
    G = M[d:N_t,:]

    return G[:,:N]*sqrt(B)
