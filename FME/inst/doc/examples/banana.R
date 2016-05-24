## =============================================================================
## A banana-shaped function.
##
## example from Haario et al. 2005
##
## Checked 24-08-09
## =============================================================================

## banana function
Banana <- function (x, a = 1, b = 1) {
  y1<-x[1]*a
  y2<-x[2]/a -b*(y1^2+a^2)
  return(c(y1,y2))
}

## banana matix function
BananaM <- function (X, a=1, b=1) {
  y1<-X[,1]*a
  y2<-X[,2]/a -b*(y1^2+a^2)
  return(cbind(y1,y2))
}


##-----------------------------------------------------------------
## Probability of a multinormally distributed value
##-----------------------------------------------------------------

pmultinorm <- function(value,mean,Cov) {
  diff <- value - mean
  ex   <- -0.5*t(diff) %*% solve(Cov) %*% diff
  rdet   <- sqrt(det(Cov))
  power  <- -length(diff)*0.5
  return((2.*pi)^power / rdet * exp(ex))
}


## 2-*log sum of squares
BananaSS <- function (p) {
  P <- Banana(p)
  Cov <- matrix(nr=2,data=c(1,0.9,0.9,1))
  -2*sum(log(pmultinorm(P,mean=0,Cov=Cov)))
}


## The Markov chain - simple Metropolis
MCMC <- modMCMC(f=BananaSS, p=c(0,0.5), jump=diag(nrow=2,x=5),niter=1000)
par(mfrow=c(4,2))
plot(MCMC,mfrow=NULL,main="MH")

## The adaptive Metropolis - update proposal every 100 runs
MCMC2 <- modMCMC(f=BananaSS, p=c(0,0.5), jump=diag(nrow=2,x=5),
                 updatecov=100,niter=1000)
plot(MCMC2,mfrow=NULL,main="AM")

## The Metropolis with delayed rejection
MCMC3 <- modMCMC(f=BananaSS, p=c(0,0.5), jump=diag(nrow=2,x=5),
                 ntrydr=2,niter=1000)

plot(MCMC3,mfrow=NULL,main="DR")

## The adaptive Metropolis with delayed rejection
MCMC4 <- modMCMC(f=BananaSS, p=c(0,0.5), jump=diag(nrow=2,x=5),
                 updatecov=100,ntrydr=2,niter=1000)

plot(MCMC4,mfrow=NULL,main="DRAM")

par(mfrow=c(2,2))
plot(MCMC$pars,main="MH")
plot(MCMC2$pars,main="AM")
plot(MCMC3$pars,main="DR")
plot(MCMC4$pars,main="DRAM")

## summaries:
colMeans(BananaM(MCMC4$pars))    # was:c(0,0)
sd(BananaM(MCMC4$pars))          # was:1
cor(BananaM(MCMC4$pars))         # 0.9 off-diagonal
