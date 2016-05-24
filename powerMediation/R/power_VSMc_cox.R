# implmement sample size calculation methods proposed in
#  Vittinghoff E., Sen, S., and McCulloch, C.E. (2009). Sample size calculation for evaluating mediation. Statistics in Medicine. 28:541-557

# for cox regression model
# b2 is the regression coefficient of mediator
# sigma.m is the standard deviation of mediator
# corr.xm is the correlation between predictor and mediator
# psi is the probability that an observation is uncensored, so that
#  the number of event d= n * psi.
ssMediation.VSMc.cox<-function(power, b2,
  sigma.m, psi, corr.xm, n.lower=1, n.upper=1e+30,
  alpha = 0.05, verbose=TRUE)
{
  if(n.lower< 1)
  {
    stop("n.lower must be >= 1")
  }
  if(n.lower >= n.upper)
  {
    stop("n.lower must be < n.upper")
  }
  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }
  if(corr.xm < -1  || corr.xm > 1)
  {
    stop("corr.xm must be in the range [-1, 1]")
  }

  res.uniroot<-uniroot(f=tmpSS.mediation.VSMc.cox,
      interval=c(n.lower, n.upper),
      power=power, b2=b2,
      sigma.m=sigma.m,
      psi=psi, corr.xm=corr.xm, alpha=alpha, verbose=FALSE)

  n.numeric<-res.uniroot$root

  res<-list(n=n.numeric, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$n) }
  invisible(res)

}

# formula (13) in Vittinghoff, Sen, and McCulloch (2009). Statist. Med. 28:541-557
# alpha - Type I error rate for one-sided test (half of Type I error rate
#    for two-sided test)
ss.VSMc.cox<-function(power, b2,
  sigma.m, psi, corr.xm, alpha = 0.025, verbose=TRUE)
{

  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }
  if(corr.xm < -1  || corr.xm > 1)
  {
    stop("corr.xm must be in the range [-1, 1]")
  }
  za<-qnorm(1-alpha)
  zg<-qnorm(power)
  numer<-(za+zg)^2
  denom<-(b2*sigma.m)^2*psi*(1-corr.xm^2)
  n.numeric<-numer/denom

  if(verbose)
  { print(n.numeric) }
  invisible(n.numeric)

}


powerMediation.VSMc.cox<-function(n, b2, sigma.m, psi, corr.xm,
  alpha=0.05, verbose=TRUE)
{
  if(n < 1)
  {
    stop("n must be >= 1")
  }
  if(corr.xm < -1  || corr.xm > 1)
  {
    stop("corr.xm must be in the range [-1, 1]")
  }

  alpha2<-alpha/2

  za2<-qnorm(1-alpha2)
  delta<-b2*sigma.m*sqrt(n*(1-corr.xm^2))*sqrt(psi)

  power<-2-pnorm(za2-delta)-pnorm(za2+delta)

  res<-list(power=power, delta=delta)
  if(verbose)
  { print(res$power) }
  invisible(res)
}


#########################
tmpSS.mediation.VSMc.cox<-function(n, power, b2,
   sigma.m, psi, corr.xm,
   alpha=0.05, verbose=FALSE)
{
  tmppower<-powerMediation.VSMc.cox(n=n, b2,
    sigma.m=sigma.m, psi=psi, corr.xm=corr.xm,
    alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}


############################
tmp.minEffect.VSMc.cox<-function(b2, n, power, sigma.m, psi, 
   corr.xm, alpha=0.05, verbose=FALSE)
{

  tmppower<-powerMediation.VSMc.cox(n=n, b2=b2, sigma.m=sigma.m, 
    psi=psi, corr.xm=corr.xm, alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}

#####################################################
minEffect.VSMc.cox<-function(n, power, sigma.m, psi, corr.xm, 
  alpha=0.05, verbose=TRUE)
{
  if(n <= 2)
  {
    stop("n must be > 2")
  }
  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }  

  res.uniroot<-uniroot(f=tmp.minEffect.VSMc.cox, 
      interval=c(0.0001, 1.0e+30),
      n=n, power=power, sigma.m=sigma.m, 
      psi=psi, corr.xm=corr.xm, alpha=alpha, verbose=FALSE)

  b2<-res.uniroot$root

  res<-list(b2=b2, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$b2) }
  invisible(res)
}

##########################


