load.Hst.ls.2Zs <-
function(Z, Z.na, Hst.ls.Z, xwhich, rgr.lags=c(0)) {
    tau <- nrow(Z)
    min.ndx <- max(1,   -min(rgr.lags)+1)
    max.ndx <- min(tau,   tau-max(rgr.lags))
	for(i in min.ndx:max.ndx) {
		
		zi.na <- Z.na[ i, ]
		Hst.ls.Z[[i]][  !zi.na , 2*xwhich-1 ] <- Z[i+rgr.lags, !zi.na ]
		Hst.ls.Z[[i]][  zi.na , 2*xwhich-1 ] <- 0
		
		Hst.ls.Z[[i]][  zi.na , 2*xwhich ] <- Z[i+rgr.lags, zi.na ]
		Hst.ls.Z[[i]][  !zi.na , 2*xwhich ] <- 0

	}
	return(Hst.ls.Z)
}
