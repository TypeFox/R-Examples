UpdateTau <-
function(survObj, priorPara, ini){

	n				<- survObj$n
	p				<- survObj$p	
	
	lambda1Sq	<- ini$lambda1Sq
	sigmaSq		<- ini$sigmaSq
	be.ini		<- ini$beta.ini
	tauSq		<- ini$tauSq


	nu.ind<-NULL
	nu=sqrt(lambda1Sq * sigmaSq/(be.ini^2))
	nu.ind <- which(nu == Inf)
	if(length(nu.ind) > 0){nu[nu.ind] <- max(nu[-nu.ind]) + 10}

	gam <- c()

	for (j in 1:p){
  		repeat{
    			gam[j]  <- rinvGauss(1, nu = nu[j], lambda = lambda1Sq)
     			if (gam[j] > 0) break    	
     	  	}
		tauSq[j] <- 1/gam[j]
		}
		
	return(tauSq)	
	
	}

