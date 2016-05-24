# $Id: dirichlet.R 2020 2015-05-23 22:12:57Z warnes $

# Posted by Ben Bolker to R-News on Fri Dec 15 2000
# http://www.r-project.org/nocvs/mail/r-help/2000/3865.html
#
# Some code (originally contributed by Ian Wilson
# <i.wilson@maths.abdn.ac.uk>


#  functions for the "Dirichlet function", the multidimensional
#  generalization of the beta distribution: it's the Bayesian
#  canonical # distribution for the parameter estimates of a
#  multinomial distribution.

# "pdirichlet" and "qdirichlet" (distribution function and quantiles)
# would be more difficult because you'd first have to decide how to
# define the distribution function for a multivariate distribution
# ... I'm sure this could be done but I don't know how



ddirichlet<-function(x,alpha)
## probability density for the Dirichlet function, where x=vector of
## probabilities
## and (alpha-1)=vector of observed samples of each type
## ddirichlet(c(p,1-p),c(x1,x2)) == dbeta(p,x1,x2)
{

  dirichlet1 <- function(x, alpha)
    {
      logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
      s <-(alpha-1)*log(x)
      s <- ifelse(alpha==1 & x==0, -Inf, s)
      exp(sum(s)-logD)
    }

  # make sure x is a matrix
  if(!is.matrix(x))
    if(is.data.frame(x))
      x <- as.matrix(x)
    else
      x <- t(x)

  if(!is.matrix(alpha))
    alpha <- matrix( alpha, ncol=length(alpha), nrow=nrow(x), byrow=TRUE)

  if( any(dim(x) != dim(alpha)) )
    stop("Mismatch between dimensions of 'x' and 'alpha'.")

  pd <- vector(length=nrow(x))
  for(i in 1:nrow(x))
    pd[i] <- dirichlet1(x[i,],alpha[i,])

  # Enforce 0 <= x[i,j] <= 1, sum(x[i,]) = 1
  pd[ apply( x, 1, function(z) any( z <0 | z > 1)) ] <- 0
  pd[ apply( x, 1, function(z) all.equal(sum( z ),1) !=TRUE) ] <- 0
  pd
}


rdirichlet<-function(n,alpha)
## generate n random deviates from the Dirichlet function with shape
## parameters alpha
{
    l<-length(alpha);
    x<-matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE);
    sm<-x%*%rep(1,l);
    x/as.vector(sm);
}

