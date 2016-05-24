simmix <- function(theta,intercept,x) {
#
# Function simmix.  To simulate data from a mixture of regressions
# model.
#

K      <- length(theta)
n      <- nrow(as.matrix(x))
theta  <- matrix(unlist(theta),ncol=length(theta))
nr     <- nrow(theta)
sigsq  <- theta[nr-1,]
lambda <- theta[nr,]
beta   <- theta[-c(nr-1,nr),]
yhat   <- if(intercept) cbind(1,x)%*%beta else x%*%beta
state  <- sample(1:K,n,TRUE,lambda)
yhat   <- yhat[n*(state-1)+(1:n)]
errr   <- rnorm(n,0,sqrt(sigsq[state]))

list(x=x,y=drop(yhat+errr))
}
