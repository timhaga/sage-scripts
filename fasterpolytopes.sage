def pairs(N,d):
    for i in xrange(1,N-d+1):
        for j in xrange(1,d+1):
            if (i,j) not in [(N-d,1), (1,d)]:
                yield (i,j)

def cover_relations(N,d):
    for i in xrange(1,N-d+1):
        for j in xrange(1,d+1):
            if (i,j) in [(1,d), (N-d,1)]:
                continue
            if i != N-d and (i,j) != (N-d-1,1):
                yield ((i,j), (i+1,j))
            if j != d and (i,j) != (1,d-1):
                yield ((i,j), (i,j+1))

def poset(N,d):
    return Poset((list(pairs(N,d)), list(cover_relations(N,d))), cover_relations=True)

def ineqs(N,d):
    # covering relations
    ineqs = [[0] + [ZZ(j==v)-ZZ(j==u) for j in pairs(N,d)]
             for u,v in cover_relations(N,d)]
    # upper bound
    ineqs += [[N] + [-ZZ(j==(N-d,d)) for j in pairs(N,d)]]
    # lower bound
    ineqs += [[0] + [ ZZ(j==(1,1)) for j in pairs(N,d)]]
    return ineqs

def eqns(N,d):
    # sum conditions
    return [[-(d*k-max(0,k-N+d)*N)] +
            [ZZ(j-i==d-k) for (i,j) in pairs(N,d)]
            for k in xrange(2,N-1)]

def posetPolytope(N,d):
    return Polyhedron(ieqs=ineqs(N,d),
            backend='ppl', base_ring=QQ)

def ambient(N,d):
    return Polyhedron(eqns=eqns(N,d),
            backend='ppl', base_ring=QQ)

def eigenstepPolytope(N,d):
    return Polyhedron(eqns=eqns(N,d), ieqs=ineqs(N,d),
            backend='ppl', base_ring=QQ)
