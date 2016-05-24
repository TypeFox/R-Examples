EstCloudiness <- function (Tx=(-999), Tn=(-999), trans=NULL, transMin = 0.15, transMax = 0.75, opt = "linear") 
{
    suppressWarnings(if ((Tx == -999 | Tn == -999) & is.null(trans)){ 
		print("Error: Please enter either Max&Min temp or transmissivity")
		} else {
	if (is.null(trans))	trans <- transmissivity(Tx, Tn)
    if (opt=="Black") {
		cl <- (0.34 - sqrt(0.34^2 + 4*0.458*(0.803-trans)))/(-2*0.458)
		cl[which(trans > 0.803)] <- 0
	} else {
		cl <- 1 - (trans-transMin) / (transMax-transMin)
	}
	cl[which(cl > 1)] <- 1
    cl[which(cl < 0)] <- 0
    return(cl)
	} )
}

