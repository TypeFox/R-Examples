UpdateLambda1.EN <-
function(survObj, priorPara, ini){

	p		<- survObj$p
	tauSq	<- ini$tauSq 
	
	r1		<- priorPara$r1
	delta1	<- priorPara$delta1
	
	lambda1Sq	<- rgamma(1, shape = r1 + p/2, rate = delta1 + sum(tauSq)/2)

	return(lambda1Sq)

	}

