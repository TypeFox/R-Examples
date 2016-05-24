
# Bayesian sample size for single arm single stage binomial trial
bayes_binom_one_postprob_onestage=function(p0,p1,eta,zeta,prior.a=0,prior.b=0,round=TRUE){

  n=3
  success=1
  failure=max(success-1,0)
  a=0;b=0
  while(a<eta | b<zeta){
    n=n+1
    a=1-pbeta(p0,prior.a+success,prior.b+n-success)
    while(a<eta){
      success=success+1
      failure=success-1
      a=1-pbeta(p0,prior.a+success,prior.b+n-success)
    }
    a=1-pbeta(p0,prior.a+success,prior.b+n-success)
    b=pbeta(p1,prior.a+failure,prior.b+n-failure)
  }

  alpha=1-pbinom(success-1,n,p0)
  power=1-pbinom(success-1,n,p1)

  if(round==TRUE){
    eta=round(a,3)
    zeta=round(b,3)
    alpha=round(1-pbinom(success-1,n,p0),3)
    power=round(1-pbinom(success-1,n,p1),3)
  }
  return(trialDesign_binom_one(reviews=n,success=success,eta=eta,zeta=zeta,alpha=alpha,power=power,exp.p0=n,exp.p1=n,p0=p0,p1=p1))

}

# ended
