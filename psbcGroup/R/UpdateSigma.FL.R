UpdateSigma.FL <-
function(survObj, priorPara, ini){
	
	p		<- survObj$p
	be.ini	<- ini$beta.ini
	covInv	<- ini$covInv

	sh.sig     <- p/2
	rate.sig   <- 1/2*as.vector(t(be.ini)%*%covInv%*%be.ini, mode = "numeric")

	sigmaSq     <- rigamma(1, a = sh.sig, b = rate.sig)
	
	return(sigmaSq)
	
	}

