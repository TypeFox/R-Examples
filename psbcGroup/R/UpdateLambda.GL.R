UpdateLambda.GL <-
function(survObj, priorPara, ini){

	p		<- survObj$p
	K		<- priorPara$K
	tauSq	<- ini$tauSq 
	
	r		<- priorPara$r
	delta	<- priorPara$delta
	
	lambda1Sq	<- rgamma(1, shape = (p + K)/2 + r, rate = delta + sum(tauSq)/2)

	return(lambda1Sq)

	}

