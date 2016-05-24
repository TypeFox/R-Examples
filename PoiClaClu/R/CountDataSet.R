CountDataSet <- function(n,p,K,param,sdsignal){
  if(n < 4*K) stop("We require n to be at least 4*K.")  
  q0 <- rexp(p, rate = 1/25)
  isDE <- runif(p) < 0.3
  classk <- matrix(NA, nrow=K, ncol=p)
  for(k in 1:K){
    lfc <- rnorm(p, sd = sdsignal)
    classk[k,] <- ifelse(isDE, q0*exp(lfc), q0)
  }
  truesf <- runif(n)*2+.2 #size factors for training observations
  truesfte <- runif(n)*2+.2 #size factors for test observations
  conds <- sample(c(rep(1:K, 4), sample(1:K, n-4*K, replace=TRUE))) # class labels for training observations
  condste <- sample(c(rep(1:K, 4), sample(1:K, n-4*K, replace=TRUE))) # class labels for test observations
  x <- xte <- matrix(NA, nrow=n, ncol=p)
  for(i in 1:n){
    for(k in 1:K){
      if(conds[i]==k) x[i,] <- rnbinom(p, mu = truesf[i]*classk[k,], size=param)
      if(condste[i]==k) xte[i,] <- rnbinom(p, mu = truesfte[i]*classk[k,], size=param)
    }
  }
  rm <- apply(x,2,sum)==0
  return(list(x=x[,!rm],xte=xte[,!rm],y=conds,yte=condste, truesf=truesf, truesfte=truesfte))
}

#CountDataSet <- function(n,p,K,param,sdsignal){
#  if(n < 4*K) stop("We require n to be at least 4*K.")  
#  q0 <- rexp(p, rate = 1/250)/10
#  isDE <- runif(p) < 0.3
#  classk <- matrix(NA, nrow=K, ncol=p)
#  for(k in 1:K){
#    lfc <- rnorm(p, sd = sdsignal)
#    classk[k,] <- ifelse(isDE, q0*2^(lfc/2), q0)
#  }
#  truesf <- pmax(.6+.3*rnorm(n), .2) #size factors for training observations
#  truesfte <- pmax(.6+.3*rnorm(n), .2)  #size factors for test observations
#  conds <- sample(c(rep(1:K, 4), sample(1:K, n-4*K, rep=TRUE))) # class labels for training observations
#  condste <- sample(c(rep(1:K, 4), sample(1:K, n-4*K, rep=TRUE))) # class labels for test observations
#  x <- xte <- matrix(NA, nrow=n, ncol=p)
#  for(i in 1:n){
#    for(k in 1:K){
#      if(conds[i]==k) x[i,] <- rnbinom(p, mu = truesf[i]*classk[k,], size=param)
#      if(condste[i]==k) xte[i,] <- rnbinom(p, mu = truesfte[i]*classk[k,], size=param)
#    }
#  }
#  rm <- apply(x,2,sum)==0
#  return(list(x=x[,!rm],xte=xte[,!rm],y=conds,yte=condste, truesf=truesf, truesfte=truesfte))
#}
