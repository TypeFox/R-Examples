N2.cohen.kappa <- function(mrg, k1, k0, alpha=0.05, power=0.8, twosided=FALSE)
{
  #	Internal function to find Maximum of St. Error
  MaxTau <- function(Marg, Kp)
  {
  	Kappa <- Kp
  	Mrg	 <- Marg
  	Dim <- length(Mrg)
  	PIe=sum(Mrg*Mrg)
  	PI0 = Kappa + PIe*(1-Kappa)
  	if(!PI0>0) stop("Invalid set of inputs", call.=FALSE)	
  	Part1=PI0*(1-PIe)^2
  	Part3=(PI0*PIe-2*PIe+PI0)^2
  	m=matrix(0,nrow=Dim, ncol=Dim)
  	for(ii in 1:Dim){
  		for (jj in 1:Dim){
  			if (ii==jj) { m[ii,ii] = -1*(1-PI0)*2*Mrg[ii]*(2*(1-PIe)-(1-PI0)*2*Mrg[ii]) }
  			else	{ m[ii,jj] = (1-PI0)^2*(Mrg[ii]+Mrg[jj])^2	}	
  		}
  	}
  	AA <- diag(1,nrow=Dim^2, ncol=Dim^2)
  	BB <- matrix(0, Dim, Dim^2)
  	CC <- diag(1,Dim,Dim)
  	for(i in 1:Dim){
  		BB[i, ((i-1)*Dim+1):(i*Dim)]=1
  		if (i<Dim) { CC <- cbind(CC, diag(1,Dim,Dim)) }
  	}
  	DD=rep(1, Dim^2)
  	EE <- as.matrix(t(as.vector(diag(1,Dim,Dim))))
  	A <- rbind(AA, BB, CC, as.matrix(t(DD)), EE)
  	rm(AA, BB, CC, DD, EE)
  	f.con <- A
  	f.obj <- as.vector(m)
  	f.dir <- c(rep(">=", Dim^2), rep("=", 2*(Dim+1)))
  	f.rhs <- c(rep(0, Dim^2), rep(Mrg, 2),1, PI0)
  	Part2=lp ("max", f.obj, f.con, f.dir, f.rhs)$objval
    #s=lp ("max", f.obj, f.con, f.dir, f.rhs)$solution
    #ss=matrix(s, byrow=T, nrow=4)
    #lp ("max", f.obj, f.con, f.dir, f.rhs)$constraints
  	TauSq <- (Part1 + Part2 - Part3)
  	if(!TauSq>0) stop("Invalid set of inputs", call.=FALSE)
  	Tau <- sqrt(TauSq)/(1-PIe)^2
  	return(Tau)
  }
  
  # ****************************************************
  # N2.cohen.kappa function
  # check for the required package
	if(!require("lpSolve"))
  		{stop("Package `lpSolve' is required and must be installed.\n 
       See help(install.packages) or write the following command at prompt
       and then follow the instructions:\n
       > install.packages(\"lpSolve\")", call.=FALSE) }
	if(any(mrg<0)) stop("Atleast one marginal probability is negative", call.=FALSE)
	if(abs(k0)>1)  stop("Invalid Value for Kappa under H0", call.=FALSE)
	if(abs(k1)>1)  stop("Invalid Value for Kappa under H1", call.=FALSE)
	if(!abs(k1-k0)>0) stop("Kappa under H0 must be different from H1", call.=FALSE)
	if(missing(alpha))	{alpha = 0.05} else {alpha = alpha}
	if(missing(power))	{power = 0.8}  else {power = power}
  if (twosided == FALSE) { side = 1 } else { side = 2 }
	if(!length(mrg)<=10) stop("Valid only up to 10 categories", call.=FALSE)
	if(!length(mrg) > 1) stop("Valid only for greater than equal to 2 categories", call.=FALSE)
	if(!sum(mrg)==1) stop("Sum of marginal probabilities is not 1", call.=FALSE)
	if(alpha<=0 || alpha>=1) stop("Alpha must be between 0 and 1", call.=FALSE)
	if(power<=0 || power>=1) stop("Power must be between 0 and 1", call.=FALSE)
	if(!(power-alpha)>0) warning("Power is less than Alpha", call.=FALSE)
	tau_Null <- MaxTau(Marg=mrg, Kp=k0)
	tau_Alt  <- MaxTau(Marg=mrg, Kp=k1)
	zAlpha   <- qnorm(1-alpha/side)
	zBeta    <- qnorm(power)
	raw_N    <- (((zAlpha*tau_Null) + (zBeta*tau_Alt))/(k0-k1))^2  
	N 	     <- ceiling(raw_N)
#	result   <- data.frame("tau_null"=tau_Null, "tau_alt"=tau_Alt, "n"=N)
#row.names(result)=""
return(N)
}
