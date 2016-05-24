Hals.fastcv.snow <-
function(j, rm.ndx, Z, Hs, Ht, Hst.ls, GP.mx) {
    n <- ncol(Z)
	tau <- nrow(Z)
    
	rho <- GP.mx[j, 1]
	reg <- GP.mx[j, 2]
    
    Z.hat <- matrix(NA, tau, n)
    
	for( drop.ndx in rm.ndx ) {
        
        if( !is.null(Hst.ls) ) {
            red.Hst.ls <- list()
            Hst0.ls <- list()
            for(i in 1:tau) { 
                red.Hst.ls[[i]] <- Hst.ls[[i]][ -drop.ndx, , drop=FALSE]
                Hst0.ls[[i]] <- Hst.ls[[i]][ drop.ndx, , drop=FALSE] 
            }
        } else {
            red.Hst.ls <- NULL
            Hst0.ls <- NULL
        }
        
        if( !is.null(Hs) ) {
            red.Hs <- Hs[ -drop.ndx, , drop=FALSE]
            Hs0 <- Hs[ drop.ndx, , drop=FALSE]
        } else {
            red.Hs <- NULL
            Hs0 <- NULL
        }
        
        Z.hat[ ,drop.ndx] <- H.als.b(Z[ , -drop.ndx, drop=FALSE], Hs=red.Hs, Ht=Ht, Hst.ls=red.Hst.ls, rho, reg, b.lag=0, Hs0=Hs0, Ht0=Ht, Hst0.ls=Hst0.ls)$Z0.hat
    }
    
	return(Z.hat)
}
