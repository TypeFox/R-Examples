## =============================================================================
## chemical example: kinetics of the reaction A<->B with forward reaction rate
## k1, and backward reaction rate k2.
## the equations: dA/dt=-k1*A + k2*B and dB/dt = k1*A - k2*B
## with initial conditions (at t=0): A0=1, B0=0.
## the analytical solution is:
## A(t) = k2/(k1+k2) + (A0-k2)/(k1+k2))*exp(-k1+k2)*t
##
## example from Haario et al. 2005
##
## Checked 24-08-09
## =============================================================================

require(FME)

## Analytic function calculates the values of A at "times".
Reaction <- function (k, times, A0=1) {
  fac <- k[1]/(k[1]+k[2])
  A   <- fac + (A0-fac)*exp(-(k[1]+k[2])*times)
  return(data.frame(t=times,A=A))
}

## Numerical function calculates the values of A and B by integration
Reaction_num <- function (k, times) {
  derivs <- function(t,state,k) {
    with (as.list(state),{
      rate <- -k[1]*A + k[2]*B
      return(list (c(dA=rate,dB=-rate)))
    })
  }
  ini <- c(A=1,B=0)
  out <- ode(func=derivs, y=ini,times=times,parms=k)
  return(as.data.frame(out))
}

## Both give the same result:
out <- Reaction_num(c(1,1),seq(0,10,by=0.25))
plot(out$time,out$A)
lines(Reaction(c(1,1),seq(0,10,by=0.25)))

## the data
Data <- data.frame(
  times = c(2,     4,     6,     8,     10   ),
  A     = c(0.661, 0.668, 0.663, 0.682, 0.650)
)

## the residual function...
residual <- function(k) return(Data$A - Reaction(k,Data$times)$A)

## the parameter prior...
Prior <- function(p,meanp=c(2,4),sigp=c(200,200))
    return( sum(((p-meanp)/sigp)^2 ))

## First fit the model to the data; initial guess = 0.5,0.5
Fit <- modFit(p=c(k1=0.5,k2=0.5),f=residual,lower=c(0,0),upper=c(1,1))
(sF <- summary(Fit))

## the residual error used as initial model variance
mse <- sF$modVariance

## the covariance matrix of fit used as initial proposal
Cov <- sF$cov.scaled * 2.4^2/2

## The Markov chain - simple Metropolis
MCMC <- modMCMC(f=residual, p=Fit$par, jump=Cov, lower=c(0,0),
                var0=mse,  prior=Prior,wvar0=1,  niter=3000)
plot(MCMC,Full=TRUE)
pairs(MCMC)

## The adaptive Metropolis - update proposal every 100 runs
MCMC2<- modMCMC(f=residual, p=Fit$par, jump=Cov, updatecov=100, lower=c(0,0),
                var0=mse, wvar0=1,    prior=Prior,niter=5000)
plot(MCMC2,Full=TRUE)
pairs(MCMC2)
sR <- sensRange(parInput=MCMC2$par,f=Reaction,times=seq(0,10,by=0.1))
plot(summary(sR))
points(Data)


## The adaptive Metropolis with delayer recjection - update proposal every 100 runs
MCMC3<- modMCMC(f=residual, p=Fit$par, jump=Cov, updatecov=100, ntrydr=3,
                lower=c(0,0), var0=mse, wvar0=0.1, prior=Prior, niter=3000)
plot(MCMC3,Full=TRUE)
pairs(MCMC3)

cumuplot(as.mcmc(MCMC$pars))
cumuplot(as.mcmc(MCMC2$pars))
cumuplot(as.mcmc(MCMC3$pars))
