tegarchRecursion <-
function(y, omega=0.1, phi1=0.4,
  kappa1=0.2, kappastar=0.1, df=10, skew=0.6,
  lambda.initial=NULL, c.code=TRUE, verbose=FALSE, aux=NULL)
{
if(is.null(aux)){
  aux <- NULL
  aux$iN <- length(y)
  aux$signnegy <- sign(-y)
  aux$u <- rep(0,aux$iN)
}
dfpluss1 <- (df+1)
u <- aux$u
mueps <- STmean(df, skew=skew)
lambda <- u
lambdadagg <- u
if(is.null(lambda.initial)){
  lambda[1] <- omega
}else{
  lambda[1] <- lambda.initial[1]
  lambdadagg[1] <- lambda.initial[2]
}

if(c.code){
  signarg <- u
  skewterm <- u
  tmp <- .tegarchRecursion(as.integer(aux$iN), as.numeric(omega),
    as.numeric(phi1), as.numeric(kappa1), as.numeric(kappastar),
    as.numeric(df), as.numeric(skew^2), as.numeric(dfpluss1),
    as.numeric(mueps), as.numeric(y), as.numeric(aux$signnegy),
    as.numeric(signarg), as.numeric(skewterm), as.numeric(lambda),
    as.numeric(lambdadagg), as.numeric(u))
  u <- tmp$u
  lambda <- tmp$lambda
  lambdadagg <- tmp$lambdadagg
}else{
  y2 <- y^2
  fn <- function(i){
    u[i] <<- dfpluss1*(y2[i]+y[i]*mueps*exp(lambda[i]))/(df*exp(2*lambda[i]) * skew^(2*sign( y[i]+mueps*exp(lambda[i]) )) + (y[i]+mueps*exp(lambda[i]))^2) - 1
    lambdadagg[i+1] <<- phi1*lambdadagg[i] + kappa1*u[i] + kappastar*aux$signnegy[i]*(u[i]+1)
    lambda[i+1] <<- omega + lambdadagg[i+1]
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
  result <- cbind(y,sigma,stdev,lambda,lambdadagg,u,epsilon,residstd)
}else{ result <- lambda }

return(result)
}
