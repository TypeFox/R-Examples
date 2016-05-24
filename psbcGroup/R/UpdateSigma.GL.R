UpdateSigma.GL <-
function(survObj, priorPara, ini){
	
	p			<- survObj$p
	be.normSq	<- ini$be.normSq
	tauSq		<- ini$tauSq 

	sh.sig     <- p/2
	rate.sig   <- 1/2*sum(be.normSq/tauSq)

	sigmaSq     <- rigamma(1, a = sh.sig, b = rate.sig)
	
	return(sigmaSq)
	
	}

