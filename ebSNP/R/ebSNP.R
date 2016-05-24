ebSNP <-
function(dat,T1=0.5,T2=0.5,eps=1e-3,maxstep=30){
  cvg <- apply(dat,2,sum)
  N <- dim(dat)[2]
  dat <- apply(dat,2,sort,decreasing=TRUE)
  p0 <- 0.8
  a0 <- 2
  b0 <- 10
  k <- dat[1,]
  out1 <- p0*choose(cvg,k)*beta(cvg-k+a0,b0+k)/beta(a0,b0)
  out2 <- 2*(1-p0)*choose(cvg,k)*(1/2)^cvg
  delta0 <- out1/(out2+out1)
  p0 <- sum(delta0)/N
  func <- function(x) -sum(delta0*log(beta(cvg-k+exp(x[1]),exp(x[2])+k)/beta(exp(x[1]),exp(x[2]))))
  para.est <- exp(optim(c(log(a0),log(b0)),func)$par)
  a0 <- para.est[1]
  b0 <- para.est[2]
  out1 <- p0*choose(cvg,k)*beta(cvg-k+a0,b0+k)/beta(a0,b0)
  out2 <- 2*(1-p0)*choose(cvg,k)*(1/2)^cvg
  delta1 <- out1/(out2+out1)
  s <- 1
  while ((sum(delta0-delta1)^2>eps)&(s<=maxstep)){
    delta0 <- delta1
    p0 <- sum(delta0)/N
    func <- function(x) -sum(delta0*log(beta(cvg-k+exp(x[1]),exp(x[2])+k)/beta(exp(x[1]),exp(x[2]))))
    para.est <- exp(optim(c(log(a0),log(b0)),func)$par)
    a0 <- para.est[1]
    b0 <- para.est[2]
    out1 <- p0*choose(cvg,k)*beta(cvg-k+a0,b0+k)/beta(a0,b0)
    out2 <- 2*(1-p0)*choose(cvg,k)*(1/2)^cvg
    delta1 <- out1/(out2+out1)
    s <- s+1
  }
  G <- rep(NA,N)
  G[which(1-delta1<T1)] <- 0
  G[which(1-delta1>=T2)] <- 1
  return(list(pi0.hat=sum(delta1)/N,
              alpha.hat=a0,
              beta.hat=b0,
              delta=delta1,
              G=G))
}
