ComputeOptimalLogLi.GivenNuHat <- function(x,y,betainit,nuhat,max){
# Uses nlm to solve for the MLE estimates for betas, given nu-hat (used to determine deviance resids)

# x = matrix of size nxp where i=1,...,n and j=1,...,p
# y = col vector of size n, i=1,...,n
# betainit = initial vector of betas, b0_1, ..., b0_p
# nuinit = initial nu value

OptimalLogLi <- rep(0,length(y))
iterct <- rep(0,length(y))

#create vector of ones
   if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x))

#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
   newx <- cbind(onevec,x)
   xmat <- as.matrix(newx)

for(i in 1:length(y)){
#for(i in 1:1){

# Create -logL = -logf (because considering single observation) so that we take the minimum of this function (which equals the max of logL)
   minuslogf <- function(par){-((y[i] * (xmat[i,] %*% par[1:length(betainit)])) - (nuhat * log(factorial(y[i]))) - log(computez(xmat[i,],par[1:length(betainit)],nuhat,max)))}

# Determine the MLEs
   BetaEstResult <- nlm(p=betainit,f=minuslogf)
#   if (BetaEstResult$code == 1 | BetaEstResult$code ==2) {OptimalLogLi[i] <- -BetaEstResult$min}
   OptimalLogLi[i] <- -BetaEstResult$min
}
return(OptimalLogLi)
}