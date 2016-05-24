UpdateTau.GL <-
function(survObj, priorPara, ini){
	
	lambdaSq	<- ini$lambdaSq
	sigmaSq		<- ini$sigmaSq
	tauSq		<- ini$tauSq
	be.normSq	<- ini$be.normSq
	
	K			<- priorPara$K
	groupInd	<- priorPara$groupInd
	groupNo		<- priorPara$groupNo
	
	nu.ind<-NULL
	nu=sqrt(lambdaSq * sigmaSq/be.normSq)
	nu.ind <- which(nu == Inf)
	if(length(nu.ind) > 0){nu[nu.ind] <- max(nu[-nu.ind]) + 10}

	gam <- c()

	for (j in 1:K){
  		repeat{
    			gam[j]  <- rinvGauss(1, nu = nu[j], lambda = lambdaSq)
     			if (gam[j] > 0) break    	
     	  	}
		tauSq[j] <- 1/gam[j]
		}
		
	return(tauSq)	
	
	}

