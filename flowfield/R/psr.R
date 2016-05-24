psr <- function (t,y) {
  #************************************************************************
  # This function extracts the penalized spline regression skeleton
  # 
  #   Refererences:  1. D. Ruppert, M. P. Wand and R. J. Carroll, Semiparametric Regression. 
  #                     New York, NY: Cambidge University Press, 2003.
  #
  #   Input: t - time series observation times 
  #          y - time series response values
  #
  #   Output:  The data frame skeleton.  The first column is the system-
  #            atically determined component (SDC) or the response variable
  #            at each knot. The second column is the forward response 
  #            derivative at each knot.
  #
  #*************************************************************************
  u <- length(t)
  knots <- floor(u/10)           # Determine number of knots
  kvector <- rep(0,knots)        # Initialize knot vector
  space <- (t[u] - t[1])/knots   # Determine knot spacing

  # Create vector of equally spaced knots between min(t) and max(t)
  for (i in 1:knots) {
    kvector[i] <- space/2 + (i-1)*space + min(t)
  }
  
  # ************************************************************************
  #  Creates the design matrix of u x (K+2) of piecewise linear functions 
  #  The matrix is used to create the "skeleton" using semi-parametric
  #  regression
  # ************************************************************************
  
  x <- matrix(data=0,nrow=u,ncol=knots+2)
  
  x[,1] <- 1
  x[,2] <- t
  for (i in 1:u) {
    for (j in 1:knots+2){
      if(t[i] > kvector[j-2]) {
        x[i,j] <- t[i] - kvector[j-2]
      }
    }
  }
  
  d <- diag(knots+2)
  d[1,1] <- 0
  d[2,2] <- 0
x <<- x
d <<- d
  # Determine the smoothing parameter
  lambda <- smoothp(t,y,x,d)
  
  b <- solve(t(x)%*%x+lambda^2*d,t(x))%*%y
  yhat <- x%*%b # fitted values
  
  # *************************************************************************
  # Extract the "skeleton" from the semi-parametric regression
  #
  # -> The skeleton consists of the smoothed response value at each knot
  #    and the forward response derivative. (The skeleton is stored in 
  #    skeleton)
  #
  # *************************************************************************
  
  sd <- rep(0,knots)    # Vector to hold the systematic component
  delta <-rep(0,knots)  # Vector to hold the slope going forward
  i <- 1
  
  while(i<=knots+1) {
    if(i==1){
      sd[i] <- b[i]+b[i+1]*min(t)
      delta[i] <- b[i+1]
    }
    else{
      sd[i] <- sd[i-1]+delta[i-1]*space
      delta[i] <- delta[i-1]+b[i+1]
    }
    i <- i+1
  }
  
  s <- x%*%solve(t(x)%*%x+lambda^2*d,t(x))
  dfres <- u - 2*sum(diag(s)) + sum(diag((s%%t(s))))
  
  sigma2 <- sum((y-yhat)^2)/dfres
  
  # Append the skelton with other information to be used later
  sd <- c(sd,knots,kvector[1],sigma2)
  delta <- c(delta,space,kvector[knots],u)
  
  skeleton <- data.frame(sd,delta)
  return(skeleton)
}




