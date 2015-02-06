@cached_function
def PolytopeOfEigensteps(N, d):
    # Create matrix of eigensteps
    E = matrix(SR, d, N+1)

    # Fill right triangle with N/d
    for k in xrange(N-d+1, N+1):
        for j in xrange(N-k , d):
            E[j,k] = N/d

    # Fill middle
    for k in xrange(1, N):
        s = max(0,d-k)
        e = min(d-1, N-k-1)
        for j in xrange(s, e+1):
            # Create variables for all but one entryy
            if j != e:
                E[j,k] = var("x_%d_%d" % (d-j, k))
            # Choose remaining entry to satisfy column sum conditions
            else:
                E[j,k] = k
                for i in xrange(0, d):
                    if i != j:
                        E[j,k] -= E[i,k]

    # Collect the interlacing conditions
    ineqs = Set()
    for k in xrange(0,N):
        for j in xrange(0,d):
            ineqs += Set([E[j,k+1] - E[j,k]])
        for j in xrange(0,d-1):
            ineqs += Set([E[j+1,k] - E[j,k+1]])

    # Transform the inequalities into a list as needed for Polyhedron()
    ineqs_l = []
    for ineq in ineqs:
        ineq_l = [ineq.polynomial(SR).constant_coefficient()]
        for k in xrange(1, N):
            s = max(0,d-k)
            e = min(d-1, N-k-1)
            for j in xrange(s, e):
                x = var("x_%d_%d" % (d-j, k))
                ineq_l.append(ineq.coefficient(x,1))
        ineqs_l.append(ineq_l)

    # Create the polytope
    return Polyhedron(ieqs=ineqs_l, backend='ppl', base_ring=QQ)

def anInteriorPoint(N, d):
        p = []
        for k in xrange(1, N):
            s = max(0,d-k)
            e = min(d-1, N-k-1)
            for j in xrange(s, e):
                i = d-j
                n = k
                p.append((d+n-2*i+1)/d)
        return p
