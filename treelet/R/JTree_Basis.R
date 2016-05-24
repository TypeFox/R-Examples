JTree_Basis <-
function(Zpos, T, PCidx, maxlev, all_nodes,whichsave){
    J = dim(Zpos)[1]
    m = dim(all_nodes)[2]
    nodes = all_nodes[maxlev, ]
    nodes = nodes[which(nodes > 0)]
    tmpfilts = diag(rep(1, m))
    ind = list()
    sums = matrix(rep(0, m * maxlev), ncol = m)
    difs = matrix(rep(0, m * maxlev), ncol = m)
    basis = list()
    for (lev in 1:maxlev) {
    	s = tmpfilts[Zpos[lev, ], ]
        R = T[[lev]]
	  y = t(R)%*%s
        tmpfilts[Zpos[lev, ], ] = y
        y = y[PCidx[lev, ], ]
        sums[lev, ] = y[1, ]
        difs[lev, ] = y[2, ]
        if (lev %in% whichsave){
        basis[[lev]]=t(tmpfilts)}
        else basis[[lev]]=NULL
    }

	return(list(basis=basis))
}
