`d2_Sigmastar_logtheta.exp` <-
function(logtheta,dist) {
		d2        = exp(- exp(logtheta) * dist) * (- exp(logtheta) * dist) *( 1 - exp(logtheta) * dist)
		diag(d2) = 0
		return(d2)
}

