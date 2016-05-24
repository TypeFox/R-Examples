UpdateSigma.EN <-
function(survObj, priorPara, ini){
	
	p		<- survObj$p
	be.ini	<- ini$beta.ini
	tauSq	<- ini$tauSq 
	lambda2	<- ini$lambda2

	sh.sig     <- p/2
	rate.sig   <- sum(1/2*be.ini^2*(1/tauSq+lambda2))

	sig.sq     <- rigamma(1, a = sh.sig, b = rate.sig)
	
	return(sig.sq)
	
	}

