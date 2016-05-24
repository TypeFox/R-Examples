H.als.b <-
function(Z, Hs, Ht, Hst.ls, rho, reg, b.lag=-1, Hs0=NULL, Ht0=NULL, Hst0.ls=NULL) {

	tau <- nrow(Z)
	n <- ncol(Z)
    Z.als <- matrix(NA, tau, n)

	if( is.null( Hs ) ) { use.Hs <- FALSE ; d.s <- 0 } else { use.Hs <- TRUE ; d.s <- ncol( Hs ) }
    if( is.null( Ht ) ) { use.Ht <- FALSE ; d.t <- 0 } else { use.Ht <- TRUE ; d.t <- ncol( Ht ) }
    if( is.null( Hst.ls ) ) { use.Hst.ls <- FALSE ; d.st <- 0 } else { use.Hst.ls <- TRUE ; d.st <- ncol( Hst.ls[[1]] ) }
    d <- d.s + d.t + d.st
	B <- matrix(0, tau, d)
	reg.mx <- diag( reg, d )
	LHH <- 0
	LHy <- 0
	g <- 1
	use.H0 <- !is.null(Hs0) | !is.null(Ht0) | !is.null(Hst0.ls)
    
	
	if( use.H0 ) { 
        if( !is.null(Hs0) ) { n0 <- nrow(Hs0) ; Z.als.0 <- matrix(NA, tau, n0) } else { n0 <- nrow(Hst0.ls[[1]]) }
    }
	
    
    low.ndx <- max( 1, 1-b.lag )
	top.ndx <- min( tau, tau-b.lag )
	for(i in low.ndx:top.ndx) {
		
        if( use.Ht ) { this.Ht.mx <- matrix( Ht[ i, ], n, d.t, byrow=TRUE ) } else { this.Ht.mx <- NULL }
        this.H <- cbind( Hs, this.Ht.mx, Hst.ls[[i]] )
        
  
		if( use.H0 ) { 
            if( use.Ht ) { this.Ht0.mx <- matrix( Ht0[ i, ], n0, d.t, byrow=TRUE ) } else { this.Ht0.mx <- NULL }
            this.H0 <- cbind( Hs0, this.Ht0.mx, Hst0.ls[[i]] )
        }

		this.HH <- crossprod(this.H)
        this.Hy <- crossprod(this.H, Z[ i, ])
		
		LHH <- LHH + g * ( this.HH - LHH )
		LHy <- LHy + g * ( this.Hy - LHy )
		inv.LHH <- try( solve( LHH + reg.mx ), silent=TRUE )
		if( class(inv.LHH) != "try=error" ) { inv.LHH <- inv.LHH } else { inv.LHH <- matrix(0, d, d) }
		
		B[ i, ] <- inv.LHH %*% LHy
		Z.als[ i, ] <- B[ i + b.lag, ] %*% t(this.H)
		if( use.H0 ) { Z.als.0[ i, ] <- B[ i + b.lag, ] %*% t(this.H0) } else { Z.als.0 <- NULL }
		g <- (g+rho) / (g+rho+1)
	}


	return( list("Z.hat"=Z.als, "B"=B, "Z0.hat"=Z.als.0,   "inv.LHH"=inv.LHH, "ALS.g"=g) )
}
