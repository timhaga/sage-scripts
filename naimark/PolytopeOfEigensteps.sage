@cached_function
def PolytopeOfEigensteps(N, d, norm_squared=1):
    # Matrix of variables
    x = matrix(SR, d, N+1)
    for n in xrange(0, N+1):
        for j in xrange(0, d):
            x[j,n] = var("x_%d_%d" % (d-j, n))

    # Collect the interlacing conditions
    ineqs = Set()
    for n in xrange(0,N):
        for j in xrange(0,d):
            ineqs += Set([x[j,n+1] - x[j,n]])
        for j in xrange(0,d-1):
            ineqs += Set([x[j+1,n] - x[j,n+1]])

    # Collect equalities
    eqs = Set()
    # First column is zero, last column is N/d*norm_squared
    for j in xrange(0,d):
        eqs += Set([x[j,0], x[j,N]-N/d*norm_squared])
    # The other columns sum to n*norm_squared
    for n in xrange(1,N):
        S = SR(0)
        for j in xrange(0,d):
            S += x[j,n]
        eqs += Set([S-n*norm_squared])

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

def anInteriorPoint(N, d, norm_squared=1):
        p = matrix(SR, d, N+1)
        for k in xrange(1, N+1):
            s = max(0,d-k)
            e = min(d, N-k)
            for j in xrange(s, e):
                i = d-j
                n = k
                p[j,k] = ((d+n-2*i+1)/d)*norm_squared
            for j in xrange(e, d):
                p[j,k] = N/d*norm_squared
        return p
