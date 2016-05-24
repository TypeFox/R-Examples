load.Hst.ls.Z <-
function(Z, Hst.ls.Z, xwhich, rgr.lags=c(0)) {
    tau <- nrow(Z)
    min.ndx <- max(1,   -min(rgr.lags)+1)
    max.ndx <- min(tau,   tau-max(rgr.lags))
	for(i in min.ndx:max.ndx) {
        Hst.ls.Z[[i]][ , xwhich ] <- t( Z[i+rgr.lags, ] )
	}
	return(Hst.ls.Z)
}
