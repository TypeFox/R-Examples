# implmement sample size calculation methods proposed in
#  Vittinghoff E., Sen, S., and McCulloch, C.E. (2009). Sample size calculation for evaluating mediation. Statistics in Medicine. 28:541-557



# sample size calculation based on
# Vittinghoff, Sen and McCulloch's (2009) method
# "Sample size calculations for evaluating mediation"
# Statistics in Medicine 2009 vol.28: 541-557
#
# (1) y=a0+a1*x+epsilon
# (2) y=b0+b1*x+b2*m+e, e~N(0, sigma.e^2)
#
# if sigma.e not available, we can use sigma.y to replace it
#  since based on the model (2), var(y|x, m)=var(e)


# for linear regression model
# b2 is the regression coefficient of mediator
# sigma.m is the standard deviation of mediator
# corr.xm is the correlation between predictor and mediator
# sigma.e is the marginal mean of the outcome
ssMediation.VSMc<-function(power, b2,
  sigma.m, sigma.e, corr.xm, n.lower=1, n.upper=1e+30,
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

  res.uniroot<-uniroot(f=tmpSS.mediation.VSMc,
      interval=c(n.lower, n.upper),
      power=power, b2=b2,
      sigma.m=sigma.m,
      sigma.e=sigma.e, corr.xm=corr.xm, alpha=alpha, verbose=FALSE)

  n.numeric<-res.uniroot$root

  res<-list(n=n.numeric, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$n) }
  invisible(res)

}

# formula (2) in Vittinghoff, Sen, and McCulloch (2009). Statist. Med. 28:541-557
# alpha - Type I error rate for one-sided test (half of Type I error rate
#    for two-sided test)
ss.VSMc<-function(power, b2,
  sigma.m, sigma.e, corr.xm, alpha = 0.025, verbose=TRUE)
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

  numer<-(za+zg)^2*sigma.e^2
  denom<-(b2*sigma.m)^2*(1-corr.xm^2)

  n.numeric<-numer/denom

  if(verbose)
  { print(n.numeric) }
  invisible(n.numeric)

}


powerMediation.VSMc<-function(n, b2, sigma.m, sigma.e, corr.xm,
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
 
  delta<-b2*sigma.m*sqrt(n*(1-corr.xm^2))/sigma.e

  power<-2-pnorm(za2-delta)-pnorm(za2+delta)

  res<-list(power=power, delta=delta)
  if(verbose)
  { print(res$power) }
  invisible(res)
}


#########################
tmpSS.mediation.VSMc<-function(n, power, b2,
   sigma.m, sigma.e, corr.xm,
   alpha=0.05, verbose=FALSE)
{
  tmppower<-powerMediation.VSMc(n=n, b2,
    sigma.m=sigma.m, sigma.e=sigma.e, corr.xm=corr.xm,
    alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}


############################
tmp.minEffect.VSMc<-function(b2, n, power, sigma.m, sigma.e, 
   corr.xm, alpha=0.05, verbose=FALSE)
{

  tmppower<-powerMediation.VSMc(n=n, b2=b2, sigma.m=sigma.m, 
    sigma.e=sigma.e, corr.xm=corr.xm, alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}

#####################################################
minEffect.VSMc<-function(n, power, sigma.m, sigma.e, corr.xm, 
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

  res.uniroot<-uniroot(f=tmp.minEffect.VSMc, 
      interval=c(0.0001, 1.0e+30),
      n=n, power=power, sigma.m=sigma.m, 
      sigma.e=sigma.e, corr.xm=corr.xm, alpha=alpha, verbose=FALSE)

  b2<-res.uniroot$root

  res<-list(b2=b2, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$b2) }
  invisible(res)
}

##########################


