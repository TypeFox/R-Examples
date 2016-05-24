
####################################
#### EBprior replicated generating
####################################
repYeb <- function(N.sim, loc, L, X=NULL, rho.family="rhoPowerExp", Y.family="Poisson", res.m=NULL, est="mode", beta=NULL, sigma=NULL, phi=NULL, k=1){
  n <- nrow(loc)
  U <- loc2U(loc)
  D <- cbind( rep(1, n), X )
  
  if(is.null(res.m)){
    if( any(c(is.null(beta), is.null(sigma), is.null(phi))) ){
      stop("res.m and the parameters are missing!")
    } else{
      m.avg <- beta; s.avg <- sigma; a.avg <- phi; k.avg <- k
    }
  } else {
    if(est=="mode"){
      if(is.vector(res.m$m)){
          m.avg <- findMode(res.m$m)
        } else{
          m.avg <- apply(res.m$m, 1, findMode)
          }
      s.avg <- findMode(res.m$s)
      a.avg <- findMode(res.m$a)
      k.avg <- ifelse(is.null(res.m$k), k, findMode(res.m$k))
    } else {
        if(is.vector(res.m$m)){
          m.avg <- mean(res.m$m)
        } else{
          m.avg <- rowMeans(res.m$m)
          }
      s.avg <- mean(res.m$s)
      a.avg <- mean(res.m$a)
      k.avg <- ifelse(is.null(res.m$k), k, mean(res.m$k))
    }
  }
  
  intcpt <- D%*%m.avg
  Z <- U2Z(U, c(s.avg, a.avg, k.avg), rho.family)
  A <- chol(Z)
  # A%*%rnorm(n) + intcpt
  S.mtx <- A%*%matrix(rnorm(n*N.sim),n,) + matrix(rep(intcpt,N.sim),n,)
    
  if(Y.family=="Poisson"){
    lamda.mtx <- exp(S.mtx)*matrix(rep(L, times=N.sim) ,n,)
    Y.mtx <- matrix(rpois(rep(1,length(S.mtx)),as.vector(lamda.mtx)), n,)
  } else if(Y.family=="Binomial"){
          p.mtx <- exp(S.mtx)/(1+exp(S.mtx))
          Y.mtx <- matrix(rbinom(rep(1,length(S.mtx)), rep(L, times=N.sim), as.vector(p.mtx)), n,)
          }
  Y.mtx
  
}
####################################
#### Post replicated generating
####################################
repYpost <- function(res.m, L, Y.family="Poisson"){
  S.mtx <- res.m$S
  n <- dim(S.mtx)[1]
  N.sim <- dim(S.mtx)[2]
  
  if(Y.family=="Poisson"){
    lamda.mtx <- exp(S.mtx)*matrix(rep(L, times=N.sim) ,n,)
    Y.mtx <- matrix(rpois(rep(1,length(S.mtx)), as.vector(lamda.mtx)), n,)
  } else if(Y.family=="Binomial"){
          p.mtx <- exp(S.mtx)/(1+exp(S.mtx))
          Y.mtx <- matrix(rbinom(rep(1,length(S.mtx)), rep(L, times=N.sim), as.vector(p.mtx)), n,)
          }
  Y.mtx
}
####################################
#### p/RPS calculation and plot
####################################
pRPS <- function(T.obs,T.rep)
  {
    p <- sum(T.rep>=T.obs)/length(T.rep)
    dd <- density(T.rep)
    h.max <- max(dd$y)
    if( T.obs<min(dd$x) || T.obs>max(dd$x) ) {
      h.obs <- 0
    } else {
      id <- sum(dd$x<T.obs)
      h.obs <- (dd$y[id+1]-dd$y[id])*(T.obs-dd$x[id])/(dd$x[id+1]-dd$x[id])+dd$y[id]
    }
    Rps <- h.obs/h.max
    res <- c(p,Rps)
    names(res) <- c("p-value","RPS")
    res
  }
###
plot_pRPS <- function(T.obs, T.rep, nm="x"){
ht.rep <- density(T.rep)
ht.obs <- ht.rep$y[sum(ht.rep$x <= T.obs )]
ind.max <- which(ht.rep$y == max(ht.rep$y))
plot(ht.rep$x, ht.rep$y, type="l", ylab=paste("h(",nm,")",sep=""),
    xlab=nm, main="Observed vs. Reference")
segments(T.obs,0, T.obs, ht.obs, col=2)
segments(ht.rep$x[ind.max],0, ht.rep$x[ind.max], max(ht.rep$y), col=3, lty=2)
legend("topright",c(paste(nm,".obs",sep=""),paste("max{h(",nm,")}",sep="")), col=2:3, lty=1:2)
}
####################################
#### Bayesian model checking
####################################
BMCT <- function(Y.obs, Y.rep, funcT, ifplot=FALSE)
{
  T.obs <- funcT(as.vector(Y.obs))
  T.rep <- apply(Y.rep,2,funcT)
  if(ifplot) plot_pRPS(T.obs,T.rep, "t")
  pRPS(T.obs,T.rep)
}
####################################
#### END
####################################