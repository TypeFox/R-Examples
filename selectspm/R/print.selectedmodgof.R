print.selectedmodgof<-function(x,...){
   
	# length of the sigma's vector
        lsig<- ((length(x$gof.u)-2)/2)
	
	hpps<-1:lsig
	hpcs<-(lsig+1):(2*lsig)
	
	cat("GoF (u) of the best model for each type: \n")
	print(c(x$gof.u[which.min(x$gof.u[hpps])],
	            x$gof.u[lsig+which.min(x$gof.u[hpcs])],
		    x$gof.u[-(1:(2*lsig))]))
	cat("\n")
	
	
	# names of the best gof
	name.bgof<- names(x$best.gof)
	#which is the best gof
	wbgof<- which.min(x$gof.u)
	
    if(wbgof <= lsig){
	 cat(paste("Best model: Inhomogeneous Poisson process (HPP) with bandwidth (sigma) =",x$best.sigma,"\n"))
       }
     if(wbgof > lsig & wbgof<= 2*lsig){
         cat(paste("Best model: Inhomogeneous Poisson cluster process (HPC) with bandwidth (sigma) =",x$best.sigma,"\n"))
	 cat(paste("and PC parameters rho =",x$best.model$rho, "and sigma^2 =",x$best.model$sigma2,"\n"))
	}
     if(wbgof > 2*lsig  & wbgof< (2*lsig)+2) cat(paste("Best model: Poisson cluster process (PC) with parameters rho =",x$best.model$rho, "and sigma^2 =",x$best.model$sigma2,"\n"))
     if(wbgof > (2*lsig)+1) cat("Best model: Poisson  process (P) \n")
   }