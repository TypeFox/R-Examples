`d2_Sigmastar_logb.exp` <-
function(logb,d) {
		d2 = diag(exp(logb),d)
		return(d2)
}

