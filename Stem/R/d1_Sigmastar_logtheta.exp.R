`d1_Sigmastar_logtheta.exp` <-
function(logtheta,dist) {
		d1       = exp(- exp(logtheta) * dist) * (- exp(logtheta) * dist)
		diag(d1) = 0
		return(d1)
}

