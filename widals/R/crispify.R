crispify <-
function( locs1, locs2, Z.delta, z.lags.vec, geodesic, alpha, flatten, self.refs, lags, stnd.d=FALSE, log10cutoff=-16 ) {	
    #	dyn.load("~/Files/Creations/C/widals.so")
	n1 <- nrow(locs1) #### length n
	n2 <- nrow(locs2) #### length n*k
    tau <- nrow(Z.delta)
	n.Zd <- ncol(Z.delta)
	Z.out <- rep(0, tau * n1)
	z.rep.in <- rep( 0:(n.Zd-1), length.out=n2 )
	t.start <- max(0, -min(lags))
	t.stop <- min(tau, tau-max(lags))
	t.s.s <- c(t.start, t.stop)
	if( geodesic ) { 
		rlocs1 <- pi*locs1[ , 1:2 ]/180
		rlocs2 <- pi*locs2[ , 1:2 ]/180
		Z.out <- .C("crispify", as.double(cos(rlocs1[,1])*cos(rlocs1[,2])), as.double(cos(rlocs1[,1])*sin(rlocs1[,2])), as.double(sin(rlocs1[,1])), as.double(locs1[ ,3]),
					as.double(cos(rlocs2[,1])*cos(rlocs2[,2])), as.double(cos(rlocs2[,1])*sin(rlocs2[,2])), as.double(sin(rlocs2[,1])), as.double(locs2[ ,3]),
					as.double( as.vector(Z.delta) ), as.double( Z.out ),
					as.double( alpha ), as.double( flatten ),
					as.integer( self.refs ), as.integer( length(self.refs) ), as.integer( z.lags.vec ), as.integer( z.rep.in ), as.integer(n.Zd),
					as.integer(n1), as.integer(n2), as.integer(tau), as.integer(stnd.d), as.integer(t.s.s), as.integer(geodesic), as.double(log10cutoff) )[[10]]
	} else { 
		Z.out <- .C("crispify", as.double( locs1[ , 1] ), as.double( locs1[ , 2] ), as.double( 0 ), as.double( locs1[ , 3] ),
					as.double( locs2[ , 1] ), as.double( locs2[ , 2] ), as.double( 0 ), as.double( locs2[ , 3] ),
					as.double( as.vector(Z.delta) ), as.double( Z.out ),
					as.double( alpha ), as.double( flatten ),
					as.integer( self.refs ), as.integer( length(self.refs) ), as.integer( z.lags.vec ), as.integer( z.rep.in ), as.integer(n.Zd),
					as.integer(n1), as.integer(n2), as.integer(tau), as.integer(stnd.d), as.integer(t.s.s), as.integer(geodesic), as.double(log10cutoff) )[[10]]
		
	}
    
    dim(Z.out) <- c(tau, n1)
    return(Z.out)
}
