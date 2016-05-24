beta_linear <- function(x, y, k, my, warm, lambda){

  nobs <- nrow(x)
  np <- ncol(x)

  beta <- matrix(data = 0.0, nrow = np, ncol = (k-1L))

  for( q in 1L:(k-1L) ) {
    temp <- numeric(np)
    for( ii in 1L:nobs ) {  
      for( jj in 1L:k ) {
        if( y[ii] == jj ) {
          temp <- temp + warm[ii,jj] * my[jj,q] * x[ii,]
        } else {
          temp <- temp - warm[ii,jj] * my[jj,q] * x[ii,]
        }
      }
    }

    beta[,q] <- temp / as.double(nobs) / lambda
  }

  return(list("beta" = beta[-1L,],
              "beta0" = beta[1L,]))

}

beta_kernel <- function(x, y, k, my, warm, lambda){

  nobs <- nrow(x)
  dnobs <- as.double(nobs)
  np <- ncol(x)

  beta <- matrix(data = 0.0, nrow = nobs, ncol = (k-1L))
  beta0 <- matrix(data = 0.0, nrow = 1L, ncol = (k-1L))

  for( q in 1L:(k-1L) ) {
    temp <- numeric(nrow(x))
    temp0 <- 0.0
    for( ii in 1L:nobs ) {
      for( jj in 1L:k ) {
        t1 <- warm[ii,jj] * my[jj,q]
        t2 <- (2.0*{y[ii] == jj} - 1.0)

        temp[ii] <- temp[ii] + t2*t1
        temp0 <- temp0 + t2*t1
      }
    }

    beta[,q] <- temp / dnobs / lambda
    beta0[,q] <- temp0 / dnobs / lambda
  }


  return(list("beta" = beta,
              "beta0" = beta0))

}


