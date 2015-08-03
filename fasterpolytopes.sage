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

def partitionOfPoint(N,d,p):
    ps = ['bottom', 'top']+list(pairs(N,d))
    labels = [0, N] + p
    return SetPartition([[ps[i] for i in xrange(d*(N-d)) if labels[i]==val] for val in set(labels)])

def showPoint(N,d,p,colorize=True,**options):
    #pos = poset(N,d)
    pos = Poset((['bottom', 'top']+list(pairs(N,d)),
                 [('bottom', (1,1)), ((N-d,d), 'top')]+
                   list(cover_relations(N,d))), cover_relations=True)
    graph = pos.hasse_diagram()

    opts = {}
    opts['pos'] = {(i,j):(i-j,i+j) for (i,j) in pairs(N,d)}
    opts['pos']['bottom'] = (0,1)
    opts['pos']['top'] = (N-2*d,N+1)
    opts['vertex_size'] = 700

    if colorize:
        opts['partition'] = partitionOfPoint(N,d,p)

    opts.update(options)
    Gplot = graph.graphplot(**opts)

    node_list = Gplot._nodelist
    pos_dict = Gplot._pos
    labels = [0, N] + p
    Gplot._plot_components['vertex_labels'] = [text(label, pos_dict[node], rgbcolor=(0,0,0), zorder=8) for node,label in zip(node_list,labels)]
    return Gplot.plot()
