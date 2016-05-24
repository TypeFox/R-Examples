UpdateLambda1.FL <-
function(survObj, priorPara, ini){

	p		<- survObj$p
	tauSq	<- ini$tauSq 
	
	r1		<- priorPara$r1
	delta1	<- priorPara$delta1
	
	lambda1Sq	<- rgamma(1, shape = r1 + p, rate = delta1 + sum(tauSq)/2)

	return(lambda1Sq)

	}

