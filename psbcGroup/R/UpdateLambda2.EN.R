UpdateLambda2.EN <-
function(survObj, priorPara, ini){

	p		<- survObj$p
	be.ini	<- ini$beta.ini
	sigmaSq	<- ini$sigmaSq 
	
	r2		<- priorPara$r2
	delta2	<- priorPara$delta2
	
	lambda2	<- rgamma(1, shape = r2 + p/2, rate = delta2 + sum(be.ini^2)/(2*sigmaSq))

	return(lambda2)

	}

