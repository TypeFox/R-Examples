UpdateLambda2.FL <-
function(survObj, priorPara, ini){

	p		<- survObj$p
	wSq		<- ini$wSq
	
	r2		<- priorPara$r2
	delta2	<- priorPara$delta2
	
	lambda2Sq	<- rgamma(1, shape = r2 + p - 1, rate = delta2 + sum(wSq)/2)

	return(lambda2Sq)

	}

