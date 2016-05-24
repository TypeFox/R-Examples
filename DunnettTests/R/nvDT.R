nvDT <- function(ratio,power,r,k,mu,mu0,contrast,sigma=NA,dist,alpha=0.05,mcs=1e+05,testcall){
  if(contrast=="means" & is.na(sigma)){
    stop("Please assume a population variance")
  }
  corr <- ratio/(1+ratio)
  if(dist=="zdist"){
    if(testcall=="SD"){
    	cvSet <- cvSDDT(k=k, alpha=alpha, alternative="U", corr=corr)
    }
    if(testcall=="SU"){
    	cvSet <- cvSUDT(k=k, alpha=alpha, alternative="U", corr=corr)
    }
  }
  
  #function to find the sample size
	samplesize <- function(n.total){
		n <- n.total/(k+1/ratio)
		n0 <- n/ratio
		if(dist=="tdist"){
			if(testcall=="SD"){
				cvSet <- cvSDDT(k=k, alpha=alpha, alternative="U",df=n*k+n0-k-1,corr=corr)
			}
      if(testcall=="SU"){
        cvSet <- cvSUDT(k=k, alpha=alpha, alternative="U",df=n*k+n0-k-1,corr=corr)
      }
		}
    
    if(contrast=="props" & is.na(sigma)){
      pooled.mu <- (mu*n*k+mu0*n0)/(n*k+n0)
      sigma <- sqrt(pooled.mu*(1-pooled.mu))
		}
    
    delta <- (mu-mu0)/(sigma*sqrt(1/n+1/n0))
    z0 <- rnorm(mcs)
		if(dist=="zdist"){u <- 1}else{u <- sqrt(rchisq(mcs,df=n*k+n0-k-1)/(n*k+n0-k-1))}
		
    dvSet <- sapply(cvSet,function(x){(x*u+sqrt(corr)*z0)/sqrt(1-corr)})        
		list.J.Fun <- list(J2.fun, J3.fun, J4.fun, J5.fun, J6.fun, J7.fun, J8.fun, J9.fun, J10.fun, J11.fun, J12.fun, J13.fun, J14.fun, J15.fun, J16.fun)
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
  pow-power
	}
    
#now solve for sample size, which should lie between 1 and the sample size required for the conservative Bonferroni
    ssBonf <- function(n.Bonf){
       if(dist=="zdist"){cv.Bonf <- qnorm(1-alpha/k)}
       if(dist=="tdist"){cv.Bonf <- qt(1-alpha/k)}
       powBonf <- 0
       n <- n.Bonf/(k+1/ratio)
	     n0 <- n/ratio
       if(contrast=="props" & is.na(sigma)){
         pooled.mu <- (mu*n*k+mu0*n0)/(n*k+n0)
         sigma <- sqrt(pooled.mu*(1-pooled.mu))
       }
       delta <- (mu-mu0)/(sigma*sqrt(1/n+1/n0))
       z0 <- rnorm(mcs)
	   if(dist=="zdist"){u <- 1}else{u <- sqrt(rchisq(mcs,df=n.Bonf-k-1)/(n*k+n0-k-1))}
       for(s in 0:(k-r)){
       	dv.Bonf <- (cv.Bonf*u+sqrt(corr)*z0)/sqrt(1-corr)
       	temp <- pnorm(dv.Bonf-delta/sqrt(1-corr))
       	powBonf <- powBonf + mean((temp^s)*((1-temp)^(k-s)))
       }
       powBonf-power
    }
  n.Bonf <- uniroot(ssBonf,interval=c(k,10000*k),tol=1e-03)$root
  n.total <- uniroot(samplesize,interval=c(k,n.Bonf),tol=1e-04)$root
  n <- ceiling(n.total/(k+1/ratio))
  n0 <- ceiling(n/ratio)
  results <- list(n,n0)
	names(results) <- c("least sample size required in each treatment groups","least sample size required in the control group")   
	return(results)
}
