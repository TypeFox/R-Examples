# Gelman-Rubin (1992) criterion for convergence (Rhat)
.GelmanRubin <- function(x,m){

  iter <- ncol(x)
  mod <- iter %% m
  n <- rep( (iter-mod)/m , m )
  nmat <- matrix(c(cumsum(n)-n+1, cumsum(n)), nrow=m)
  n <- n[1]

  Rhat <- numeric(nrow(x))
  for(ii in 1:nrow(x)){

    # values per chain
    chs <- apply(nmat, 1, function(j) x[ii,j[1]:j[2]])
    mns <- apply(chs,2,mean)
    vrs <- apply(chs,2,var)
    Bdivn <- var(mns)
    W <- mean(vrs)
    muhat <- mean(chs)
    sighat2 <- (n-1)/n * W + Bdivn
    # sampling distribution
    Vhat <- sighat2 + Bdivn/m
    var.Vhat <- ((n-1)/n)^2*(1/m)*var(vrs) + ((m+1)/(m*n))^2*2/(m-1)*(Bdivn*n)^2 +
                2*((m+1)*(n-1)/(m*n^2)) * (n/m)*(cov(vrs,mns^2)-2*muhat*cov(vrs,mns))
    df <- 2*Vhat^2 / var.Vhat
    # compute Rhat
    if(Bdivn==0 & identical(vrs,rep(0,m))){ # for zero variance defined as 1
      Rhat[ii] <- 1
    }else{
      Rhat[ii] <- sqrt( (Vhat/W)*df/(df-2) )
    }

  }
  Rhat

}

# criterion for goodness of approximation (Hoff, 2009)
.SDprop <- function(x){

  np <- nrow(x)
  v <- apply(x, 1, var)   # variance of chain
  v0 <- v==0
  sdp <- sp0 <- numeric(np)
  for(i in 1:np){
    arp <- try( ar(x[i,], aic=TRUE), silent=T )
    if(!v0[i]) sp0[i] <- arp$var.pred/(1 - sum(arp$ar))^2   # spectral density at frequency 0
  }
  mcmc.v <- sp0/ncol(x)   # mcmc-variance (correcting for autocorrelation)
  # proportion of variance due to sampling inefficiency
  sdp[!v0] <- sqrt(mcmc.v / v)[!v0] # for zero variance defined as 0
  sdp

}

