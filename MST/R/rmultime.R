rmultime <-
function(N=100, K=4, beta=c(-1, 2, 1, 0, 0), cutoff=c(.5, .5, 0, 0),
                     digits=1, icensor=1,
                     model = c("gamma.frailty", "log.normal.frailty",
                               "marginal.multivariate.exponential", "marginal.nonabsolutely.continuous",
                               "nonPH.weibull"),
                     v=1, rho=.65, a=1.5, lambda=0.1){
  model<-match.arg(model,c("gamma.frailty", "log.normal.frailty","marginal.multivariate.exponential",
                           "marginal.nonabsolutely.continuous", "nonPH.weibull"))

  #### Generate Covariates
  p <- length(beta)-1 
  X <- matrix(0, nrow=N*K, ncol=p)
  for (j in 1:p) {
    if (is.odd(j)) X[,j] <- rep(runif(N, 0, 1), rep(K,N))
    else X[, j] <- runif(N*K)
  }
  X <- round(X, digits=digits)
  Z <- X
  for (j in 1:length(cutoff)) {
    if (cutoff[j] > 0 && cutoff[j] <1) Z[,j] <- sign(X[,j] <= cutoff[j])
  }
  Z <- cbind(1, Z)
  eta <- Z%*%beta	
  
  #### Generate Multivariate Survival and Censoring Time Data from Different Models
  if (model=="gamma.frailty") {
    if (v==0){w <- 1
    } else {w <- rep(rgamma(N, 1/v, 1/v), rep(K, N))}
    rate <- exp(eta)*w 
    xobs <- rexp(N*K, rate)
    dind <- 1
    if (icensor!=0){     
      c <- rexp(N*K, rate)
      dind <- sign(xobs <= c*icensor)
      xobs <- pmin(xobs, c*icensor) 
    } 
  } else if (model=="log.normal.frailty") {
    w <- 1
    if ( v!=0) w <- rep(rnorm(N, 0, v), rep(K, N))
    rate <- exp(eta + w)
    xobs <- rexp(N*K, rate)
    dind <- 1
    if (icensor!=0){     
      c <- rexp(N*K, rate)
      dind <- sign(xobs <= c*icensor)
      xobs <- pmin(xobs, c*icensor) 
    }
  } else if (model=="marginal.multivariate.exponential") {
    cor <- matrix(rho, K, K)
    diag(cor) <- 1
    x <- -log(1-pnorm(as.vector(t(mvrnorm(N,mu=rep(0,K),Sigma=cor)))))
    c <- -log(1-pnorm(as.vector(t(mvrnorm(N,mu=rep(0,K),Sigma=cor)))))
    if (icensor == 0) {
      xobs <- x* exp(-eta)
      dind <- 1
    } else {
      dind <- sign(x <= c*icensor)
      xobs <- pmin(x, c*icensor)*exp(-eta)
    }    
  } else if (model =="marginal.nonabsolutely.continuous"){
    lambda0 <- 2 - 2/(rho+1)    
    x <- pmin(rep(rexp(N, lambda0), rep(K,N)),  rexp(N*K, 1-lambda0))
    c <- pmin(rep(rexp(N, lambda0), rep(K,N)),  rexp(N*K, 1-lambda0))
    if (icensor==0) {
      xobs <- x* exp(-eta)
      dind <- 1
    } else {     
      xobs <- pmin(x, c*icensor)* exp(-eta)
      dind <- sign(x <= c*icensor)
    }
  } else if (model=="nonPH.weibull"){
    if (v==0) w <- 1
    else {w <- rep(rgamma(N, 1/v, 1/v), rep(K, N))}
    rate <- exp(eta)*w
    lambda1 <- lambda*rate
    b <- (lambda1)^-(1/a) 
    xobs <- rweibull(N*K, shape=a, scale=b) 
    dind <- 1
    if (icensor!=0) {
      c <- rweibull(N*K, shape=a, scale=b) 
      dind <- sign(xobs <= c*icensor)
      xobs <- pmin(xobs, c*icensor) 
    }
  }
  
  ##### Output
  colnames(X) <- paste("x", 1:p, sep="")
  dat <- data.frame(id=rep(1:N, rep(K,N)),rep=rep(1:K, N), 
                    time=xobs, status=dind, X)
  return(list(dat=dat, model=model))			  
}
