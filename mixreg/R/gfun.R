gfun <- function(x,y,theta) {
#
# Function gfun.  To calculate the probabilites gamma[i,j]
# that observation i corresponds to state j for the Gaussian model.
#

# Argument theta is a list of lists; theta[[j]] has
# components beta, sigsq, and lambda.

K      <- length(theta)
theta  <- matrix(unlist(theta),ncol=K)
nr     <- nrow(theta)
sigsq  <- theta[nr-1,]
lambda <- theta[nr,]
beta   <- theta[-c(nr-1,nr),]
n      <- nrow(x)

yhat   <- x%*%beta
fff    <- matrix(dnorm(t(y-yhat)/sqrt(sigsq)),K,n)/sqrt(sigsq)

hhh    <- t(lambda*fff)
ll     <- apply(hhh,1,sum)
list(gamma=hhh/ll,log.like=sum(log(ll)))
}
