unload.Hst.ls <-
function(Hst.ls, which.col, rgr.lags) {
	n <- nrow(Hst.ls[[1]])
	tau <- length(Hst.ls)
	Z.out <- matrix(NA, tau, n)
    min.ndx <- max(1,   -min(rgr.lags)+1)
    max.ndx <- min(tau,   tau-max(rgr.lags))
	for(i in min.ndx:max.ndx) {
		Z.out[i-rgr.lags, ] <- Hst.ls[[i]][ , which.col ]
	}
	return(Z.out)
}
