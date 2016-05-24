optim.ftn.2 <-
function(dat,coef, random.eff,iter){
  n.coef=length(coef)
  for(pp in 1:iter){     
    err=0  
    err.num=0
    result1=try(optim.step1.2(dat, coef, random.eff), silent=T)
    if(class(result1)== "try-error"){err=1}
    while(err==1){
      coef[3]=coef[3]+0.002
      result1=try(optim.step1.2(dat, coef, random.eff), silent=T)
      if(err.num>3){return(res=list(coef=rep(0, n.coef), ran.eff=NULL,fitted=NULL, yhat=NULL, 
                                    loglik=NULL, error=NULL, std.error=NULL, nlme.result=NULL, dat=dat, conv=F))
                    break}
      err.num=err.num+1
      if(class(result1)!= "try-error"){err=0}
    }
    result2=nlme(DAMAGE_Y~optim.step2.2(Dt,  A, B, r1, r2), fixed=A+B~1,random=r1+r2~1 | id, 
                 start=c(result1$coef[1], result1$coef[2]), data=result1$newdata)
    
    kk=as.vector(attr(result2$apVar,"Pars"))
    ss0=try(exp(kk[1]), silent=T)
    ttt=0
    
    if(class(ss0)== "try-error"){ttt=1} 
    while(ttt==1){
      tryA=rnorm(1, result1$coef[1], 0.5)
      tryB=rnorm(1, result1$coef[2], 0.5)
      result2=nlme(DAMAGE_Y~optim.step2.2(Dt,  A, B, r1, r2), fixed=A+B~1,random=r1+r2~1 | id, 
                   start=c(tryA, tryB), data=result1$newdata)
      kk=as.vector(attr(result2$apVar,"Pars"))    
      ss0=try(exp(kk[1]), silent=T)
      if(class(ss0)== "try-error") { ttt=1 }else{ttt=0}
    }
    ss0=exp(kk[1])  
    ss1=exp(kk[2])
    rho=2/(1+exp(-kk[3]))-1
    ss=exp(kk[4])
    
    ss.vec=c(ss0, ss1, rho, ss)
    A=result2$coef$fixed[1]
    B=result2$coef$fixed[2]
    random.eff=result2$coef$random$id
    beta=result1$coef[3:(n.coef-4)] 
    theta0=ss.vec
    
    diff.A=abs(A-result1$coef[1])
    diff.B=abs(B-result1$coef[2])
    diff.beta=abs(coef[3:(n.coef-4)]-beta)
    diff.random.eff=abs(random.eff-result1$random.eff) 
    diff.theta0=abs(theta0-coef[(n.coef-3):n.coef])    
    
    diff.max=max(diff.A, diff.B, diff.random.eff, diff.beta, diff.theta0)
    #print(c("diff.max=", diff.max))
    coef=c(A, B, beta, theta0)
    dfs=dat$dfs
    names(coef)=c("G", "H", colnames(dat$dat)[4:(sum(dfs)+4)], "sigma0", "sigma1", "rho", "sigma")
  }
  loglik=-minus.log.lik.nlme(dat, coef, random.eff)
  error=result1$residual
  std.error=error/coef[length(coef)]
  res=list(coef=coef, ran.eff=random.eff,fitted=result1$yhat, yhat=result1$yhat, 
           loglik=loglik, error=error, std.error=std.error, nlme.result=result2, dat=dat, conv=T)
  return(res)    
}
