tegarchRecursion2 <-
function(y, omega=0.1, phi1=0.4,
  phi2=0.2, kappa1=0.05, kappa2=0.1, kappastar=0.02, df=10,
  skew=0.6, lambda.initial=NULL, c.code=TRUE, verbose=FALSE,
  aux=NULL)
{
if(is.null(aux)){
  aux <- NULL
  aux$iN <- length(y)
  aux$signnegy <- sign(-y)
  aux$u <- rep(0,aux$iN)
}
dfpluss1 <- (df+1)
y2 <- y^2
mueps <- STmean(df, skew=skew)
u <- aux$u
lambda <- u #lambda
lambda1dagg <- rep(0,aux$iN) #lambda1dagg
lambda2dagg <- lambda1dagg #lambda2dagg
if(is.null(lambda.initial)){
  lambda[1] <- omega
}else{
  lambda[1] <- lambda.initial[1]
  lambda1dagg[1] <- lambda.initial[2]
  lambda2dagg[1] <- lambda.initial[3]
}

if(c.code){
  signarg <- u
  skewterm <- u
  tmp <- .tegarchRecursion2(as.integer(aux$iN), as.numeric(omega),
    as.numeric(phi1), as.numeric(phi2), as.numeric(kappa1), as.numeric(kappa2),
    as.numeric(kappastar), as.numeric(df), as.numeric(skew^2),
    as.numeric(dfpluss1), as.numeric(mueps), as.numeric(y), as.numeric(y2),
    as.numeric(aux$signnegy), as.numeric(signarg), as.numeric(skewterm),
    as.numeric(lambda), as.numeric(lambda1dagg), as.numeric(lambda2dagg),
    as.numeric(u))
  u <- tmp$u
  lambda1dagg <- tmp$lambda1dagg
  lambda2dagg <- tmp$lambda2dagg
  lambda <- tmp$lambda
}else{
  fn <- function(i){
    u[i] <<- dfpluss1*(y2[i]+y[i]*mueps*exp(lambda[i]))/(df*exp(2*lambda[i]) * skew^(2*sign( y[i]+mueps*exp(lambda[i]) )) + (y[i]+mueps*exp(lambda[i]))^2) - 1
    lambda1dagg[i+1] <<- phi1*lambda1dagg[i] + kappa1*u[i]
    lambda2dagg[i+1] <<- phi2*lambda2dagg[i] + kappa2*u[i] + kappastar*aux$signnegy[i]*(u[i]+1)
    lambda[i+1] <<- omega + lambda1dagg[i+1] + lambda2dagg[i+1]
  }
  indx <- 1:I(aux$iN-1)
  tmp <- sapply(indx,fn)
}

#output:
if(verbose){
  u[aux$iN] <- NA
  sigma <- exp(lambda)
  epsilon <- y/sigma
#  sdepsilon <- sd(epsilon)
  sdepsilon <- sqrt(STvar(df, skew=skew))
  stdev <- sigma*sdepsilon
  residstd <- epsilon/sdepsilon
  result <- cbind(y,sigma,stdev,lambda,lambda1dagg,lambda2dagg,u,epsilon,residstd)
}else{ result <- lambda }

return(result)
}
