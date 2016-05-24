ComputeBetasAndNuHat <- function(x,y,betainit,nuinit,max){
# Uses nlminb to solve for the MLE estimates for betas and nu

# x = matrix of size nxp where i=1,...,n and j=1,...,p
# y = col vector of size n, i=1,...,n
# betainit = initial vector of betas, b0_1, ..., b0_p
# nuinit = initial nu value

#create vector of ones
   if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x))

#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
   newx <- cbind(onevec,x)
   xmat <- as.matrix(newx)

# Create -logL so that we take the minimum of this function (which equals the max of logL)
   minusloglike <- function(par){-sum((y * (xmat %*% par[1:length(betainit)])) - (par[length(betainit)+1] * log(factorial(y))) - log(computez(xmat,par[1:length(betainit)],par[length(betainit)+1],max)))}

# Determine the MLEs
   BetaNuEst <- nlminb(start=c(betainit,nuinit),minusloglike,lower = c(rep(-Inf,length(betainit)),0), upper = c(rep(Inf,length(betainit)),Inf))

   if (BetaNuEst$convergence != 0) {
	# Determine the MLEs
	   BetaNuEst <- optim(par=c(betainit,nuinit),fn=minusloglike,control=list(maxit=1000))#$par
   }

	return(BetaNuEst)
}

