RepLOptim <- function(parmean,parsd,fr,maxrepet,maxnoimprov,maxreplic,maxiter,gr=NULL,inphess=NULL,method="nlminb",lower=NULL,upper=NULL,rethess=FALSE,
allrep=NULL,maxeval=NULL,objbnd=Inf,parmstder=FALSE,tol=1E-8,...) 
{
   npar <- length(parmean)
    values <- NULL
    if (is.null(lower)) lower <- rep(-Inf,npar)
    if (is.null(upper)) upper <- rep(Inf,npar)
    if (is.finite(objbnd))  { if (is.null(allrep)) allrep <- 10*maxreplic }
    else allrep <- maxreplic

    bestres <- NULL	
    bestval <- Inf
    bestpar <- parmean
    initpar <- parmean
    cnt <- repcnt <- noimpcnt <- 0 
    for (i in 1:maxreplic)  {
       if (cnt > allrep | repcnt >= maxrepet | noimpcnt >= maxnoimprov) break
       value <- Inf
       while (value >= objbnd && cnt < allrep)
      {
       	 if (method == "nlminb")
	    tmpres <- nlminb(start=initpar,fr,gradient=gr,hessian=inphess,lower=lower,upper=upper,
				control=list(iter.max=maxiter,eval.max=maxeval),...)
       	 else if (method == "nlm") 
	    tmpres <- nlm(fr,p=initpar,lbound=lower,ubound=upper,iterlim=maxiter,...)
       	 else if (method == "L-BFGS-B")
	    tmpres <- optim(initpar,fr,gr=gr,method=method,lower=lower,upper=upper,
				control=list(maxit=maxiter),hessian=rethess,...)
       	 else  tmpres <-
            optim(initpar,fr,gr=gr,method=method,control=list(maxit=maxiter),lbound=lower,ubound=upper,hessian=rethess,...)
       	 if (method == "nlminb") value <- tmpres$objective
       	 else if (method == "nlm") value <- tmpres$minimum
       	 else value<- tmpres$value
         if (is.na(value)) value <- objbnd
         cnt <- cnt+1
	 u <- runif(n=npar)     # generate npar uniform random numbers
       	 initpar <- qnorm(u,mean=bestpar,sd=parsd) #  generate new parameters from a normal distribution
	 lbndind<- initpar < lower   #  identify indices of parameters that fell below their lower bounds
	 ubndind <- initpar > upper   #  identify indices of parameters that fell above their upper bounds
	 initpar[lbndind] <- lower[lbndind] + u[lbndind] * (bestpar[lbndind]-lower[lbndind]) # and correct those parameters
	 initpar[ubndind] <- upper[ubndind] - u[ubndind] * (upper[ubndind]-bestpar[ubndind])
       } 
      values <- c(values,value)
      if (is.na(value)) { noimpcnt <- noimpcnt + 1 ; repcnt <- 0 }
      else {
      	if (is.finite(bestval))  if (abs((value-bestval)/bestval) < tol) repcnt <- repcnt + 1 
      	else repcnt <- 0 
      	if (value < bestval)  {
          bestval <- value
          if (method != "nlm") bestpar <- tmpres$par
          else bestpar <- tmpres$estimate
          bestres <- tmpres
      	  noimpcnt <- 0 
       }
       else noimpcnt <- noimpcnt + 1
     }
    }

    if (method == "nlminb")  {
        iterations <- bestres$iterations
        if (method != "nlm") counts <- bestres$evaluations
        else counts <- NULL
        hess <- NULL
        egval <- NULL
        parstd <- NULL
    } 
    else  {
    	iterations <- NULL
    	counts <- bestres$counts
        if (rethess==TRUE) {
           hess <- bestres$hessian
           egval <- eigen(hess,symmetric=TRUE,only.values=TRUE)$values
       	   if (parmstder==TRUE)
              if (egval[npar] < tol) parstd <- "Not computed because the hessian is not positive definite"
              else parstd <- sqrt(diag(solve(hess)))
           else parstd <- NULL
        }
        else {
	    hess <- NULL
            egval <- NULL
            parstd <- NULL
       }
    } 
    if (!is.null(bestres))
    	return(list(par=bestpar,val=bestval,iterations=iterations,vallist=values,counts=counts,convergence=bestres$convergence,message=bestres$message,hessian=hess,hessegval=egval,stderrors=parstd))
    else
    	return(list(par=NULL,val=Inf,iterations=NULL,vallist=NULL,counts=NULL,convergence=NULL,message="RepLOptim was unable to find any valid solution",hessian=NULL,hessegval=NULL,stderrors=NULL))
}
