# Frequentist sample size for single arm single stage binomial trial
freq_binom_one_onestage=function(p0,p1,alpha,power,prior.a=0,prior.b=0,round=TRUE){
  n=0
  a=0;b=0
  while(b<power){
    n=n+1
    a=qbinom(1-alpha,n,p0)
    b=1-pbinom(a,n,p1)
  }

  alpha=1-pbinom(a,n,p0)
  eta=round(1-pbeta(p0,prior.a+a+1,prior.b+n-a-1),3)
  zeta=round(pbeta(p1,prior.a+a,prior.b+n-a),3)
  success=a+1

  if(round==TRUE){
    alpha=round(alpha,3)
    power=round(b,3)
    eta=round(eta,3)
    zeta=round(zeta,3)
  }

  return(trialDesign_binom_one(reviews=n,success=success,eta=eta,zeta=zeta,alpha=alpha,power=power,exp.p0=n,exp.p1=n,p0=p0,p1=p1))

}

# ended
