UpdateW.FL <-
function(survObj, priorPara, ini){
	
	lambda2Sq	<- ini$lambda2Sq
	sigmaSq		<- ini$sigmaSq
	be.ini		<- ini$beta.ini
	p <- dim(survObj$x)[2]
	
	nu.ind<-NULL
	diff.beta <- diff(be.ini)
	
	nu=sqrt(lambda2Sq * sigmaSq/(diff.beta^2))	
	nu.ind <- which(nu == Inf)
	
	if(length(nu.ind) == p-1){
		nu <- rep(sqrt(lambda2Sq * sigmaSq/0.00000001), p-1)
		}
	if(length(nu.ind) >0 & length(nu.ind) < p-1){
		nu[nu.ind] <- max(nu[-nu.ind]) + 10
		}	
		
	wSq.inv <- c()
	wSq		<- c()

	for (j in 1:(p-1)){
  		repeat{
    			wSq.inv[j]  <- rinvGauss(1, nu = nu[j], lambda = lambda2Sq)
     			if (wSq.inv[j] > 0) break    	
     	  	}
	wSq[j] <- 1/wSq.inv[j]
	}	
	return(wSq)	
	
	}

