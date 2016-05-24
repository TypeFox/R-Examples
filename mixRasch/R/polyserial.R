'polyserial' <- function(x, y, ml=TRUE){

	# y is categorical
	# x is continuous
	if(! is.numeric(x)){
		x <- is.numeric(x)
		warning("x was coerced to numeric.")
	}
	noMissing <- complete.cases(x,y)
    x <- x[noMissing]
    z <- scale(x) 
    
    y <- y[noMissing]
    y <- as.numeric(as.factor(y))
    k <- max(y)
    N <- length(y)
    tau <- table(y)
    tau <- qnorm(cumsum(tau[-k])/N)
    
    if(k==1){
    	return(NA)
    }
    
    r.xy <- cor(x,y)
    
    s.pop <- sd(y)*sqrt((N-1)/(N))
    
    rho <- (r.xy * s.pop)/sum(dnorm(tau)) 
    
     if(ml){
      negLogLik <- function(par){
    	rho <- par[1]
		if(rho >= 1) rho <- .99999
		if(rho <= -1) rho <- -.99999
    	tau <- c(-Inf,par[-1],Inf)
    	tauDenom <- sqrt(1-rho^2)
		rhoZ <- rho*z
    	pTauStar <- pnorm((tau[y+1] - rhoZ)/tauDenom)
		pTauStarLower <- pnorm((tau[y] - rhoZ)/tauDenom) 
		pTauDif <- pTauStar-pTauStarLower
		pTauDif <- ifelse(pTauDif == 0, .5e-10, pTauDif)
    	-sum(log(dnorm(z))+log(pTauDif))
      }
      mlEst <- optim(c(rho,tau), negLogLik)	
      rho <- mlEst$par[1]
    }
    rho	
}