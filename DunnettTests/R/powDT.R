powDT <- function(r,k,mu,mu0,n,n0,contrast,sigma=NA,df=Inf,alpha=0.05,mcs=1e+05,testcall){
	#step 1: calculate the standardized effect size
	theta <- mu-mu0
	if(contrast=="props" & is.na(sigma)){
		pooled.mu <- (mu*n*k+mu0*n0)/(n*k+n0)
		sigma <- sqrt(pooled.mu*(1-pooled.mu))
	}
	corr <- n/(n0+n)
	delta <- theta/(sigma*sqrt(1/n+1/n0))
    
  #step 2: calculate the critical values 
  if(testcall=="SD"){
    cvSet <- cvSDDT(k=k, alpha=alpha, alternative="U", df=df, corr=corr)
  }
  if(testcall=="SU"){
    cvSet <- cvSUDT(k=k, alpha=alpha, alternative="U", df=df, corr=corr)
  }
    
  #step 3: numerically approximate the power by monte carlo.
  z0 <- rnorm(mcs)
  if(df == Inf){u <- 1}else{u <- sqrt(rchisq(mcs,df=df)/df)}
  dvSet <- sapply(cvSet,function(x){(x*u+sqrt(corr)*z0)/sqrt(1-corr)})
  # with columns corresponds to c1,c2,...,ck
  list.J.Fun <- list(J2.fun,J3.fun,J4.fun,J5.fun,J6.fun,J7.fun,J8.fun,J9.fun,J10.fun,J11.fun,J12.fun,J13.fun,J14.fun,J15.fun,J16.fun)
  if(testcall=="SU"){
    	#initialize the pow as the probability of rejecting all the k nulls
    	pow <- mean((1-pnorm(dvSet[,1],mean=delta/sqrt(1-corr)))^k)
    	if(r<k){
    		#add the probabilities of rejecting k-s out of k nulls, s=1,...,k-r
            s <- 1
            p.rej <- (1-pnorm(dvSet[,s+1],mean=delta/sqrt(1-corr)))^(k-s)
    		Psi.d1 <- pnorm(dvSet[,s],mean=delta/sqrt(1-corr))
    		J1 <- Psi.d1
    		list.J <- c(list(1),list(J1))
    		list.Psi <- list(Psi.d1)
    		pow <- pow+choose(k,s)*mean(p.rej*J1)
    		if(r<(k-1)){
    			for(s in 2:(k-r)){
    	       	p.rej <- (1-pnorm(dvSet[,s+1],mean=delta/sqrt(1-corr)))^(k-s)
    	       	Psi.ds <- pnorm(dvSet[,s],mean=delta/sqrt(1-corr))
    			Js <- list.J.Fun[[s-1]](Psi.ds,list.J,list.Psi)
    			list.J <- c(list.J,list(Js))
    			list.Psi <- c(list.Psi,list(Psi.ds))
    			pow <- pow+choose(k,s)*mean(p.rej*Js)
    			}
    		}
    	}
  }

  if(testcall=="SD"){
    	powSet <- NULL #with elements correspondng to reject 1, 2,...,k in sequence
    	Psi.d1 <- 1-pnorm(dvSet[,k],mean=delta/sqrt(1-corr)) 
    	J1 <- Psi.d1
    	list.J <- c(list(1),list(J1))
    	list.Psi <- list(Psi.d1)
    	p.acc <- pnorm(dvSet[,k-1],mean=delta/sqrt(1-corr))^(k-1)
    	powSet <- c(powSet,choose(k,k-1)*mean(p.acc*J1)) #reject only 1
    	if(k>2){
    		for(s in (k-2):0){
    			if(s==0){p.acc <- 1}else{p.acc <- pnorm(dvSet[,s],mean=delta/sqrt(1-corr))^s}
    			Psi.ds <- 1-pnorm(dvSet[,s+1],mean=delta/sqrt(1-corr))
    			Js <- list.J.Fun[[k-s-1]](Psi.ds,list.J,list.Psi)
    			powSet <- c(powSet,choose(k,s)*mean(p.acc*Js)) #reject k-s
    			list.J <- c(list.J,list(Js))
    			list.Psi <- c(list.Psi,list(Psi.ds))
    		}
    	}
      pow <- sum(powSet[r:k]) 
  }
pow
}
