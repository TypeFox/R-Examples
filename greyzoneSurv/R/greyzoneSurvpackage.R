
# packages=c('stats4', 'survival', 'Hmisc', 'survAUC')
# for (package in packages)
# {
#   if(package %in% rownames(installed.packages()))
#     do.call('library', list(package))
#   
#   # if package is not installed locally, download, then load
#   else {
#     install.packages(package)
#     do.call("library", list(package))
#   }
#     
# }

genSurvData=function(n, recruitment.yrs=2, baseline.hazard=365.25*5, shape=1, censoring.rate=0, beta.continuous, beta.binary=0, x, xhigh, ran.seed)
{ 
  #### This function generates weibull survival times with a fixed censoring rate.
  #### In the R rweibull function, the scale is different than the usual lamda used for the scale parameter. 
  #### In the usual parameterization, the bigger the lamda, the smaller the mean of the survival time.
  #### Here in rweibull, the bigger the scale the bigger the mean of the survival time.


  ###### simulate exponential survival times
  lambdaT=baseline.hazard
  
  #true event time for low risk group, note that lambdaT is the mean of survival times when covariate=0
    set.seed(ran.seed)
    T0 = floor(rweibull(sum(xhigh==0), shape=shape, scale=lambdaT*exp(-beta.continuous*x[xhigh==0])))
    #true event time for high risk group
    set.seed(ran.seed)  
    T1=floor(rweibull(sum(xhigh==1), shape=shape, scale=lambdaT*exp(-beta.continuous*x[xhigh==1]-beta.binary)))
    T=c(T0, T1); T
    
  if (censoring.rate==0)
  {
    event=rep(1, n)
    time=T
    
  } else
  {
    #simulate the dates when patients enter the study
    set.seed(ran.seed)  
    ondt=sample.int(ceiling(recruitment.yrs*365.25),size=n); ondt
    
    #number of events at minimum (when there are ties there might be more events therefore fewer censored)
    n.events0=n-floor(n*censoring.rate); n.events0
    #what patients to censor? Note we need to get at least n.events0 events (it could be more because of tied survival times)
    tot=ondt+T 
    tot.ordered=sort(tot)
    #want to censor after tot.ordered[n.events0], so need to randomly choose a day as stopping-study date 
    #between tot.ordered[n.events0] and tot.ordered[n.events0+1]-1
    #set.seed(ran.seed)  
    tot.ordered[n.events0]
    tot.ordered[n.events0+1]
    num.dup=sum(tot==tot.ordered[n.events0]); num.dup#this tells the number of duplicated values at where we want to stop the study
    days.to.choose=seq(tot.ordered[n.events0],tot.ordered[n.events0+num.dup]-1, by=1)
    if (length(days.to.choose)>1) stop.dt=sample(days.to.choose, size=1) else stop.dt=tot.ordered[n.events0]  #this is study-stopping date
    stop.dt
    
    event=(T+ondt<=stop.dt)*1
    time=pmin(T+ondt, stop.dt)-ondt
    
  }
  
  return(out=list(data=data.frame(time, event, x, xhigh), censoring.rate=mean(event==0)))  
  
}


em.func=function(error=.001, max.iter=300, initial.values=NULL, 
                  y, delta, x)
{
 
  #############################################
  #  EM to get parameter estimates from joint 
  # Weilbull and logistic likelihood function
  #############################################
  
  # Meanings of the function arguments:
  # z is a vector of whether a patient is high risk (=1) or low risk (=0)
  # lamda is a value
  # beta=c(beta0, beta1)
  # gamma=c(gamma0, gamma1)
  # y is a vector of survival time
  # delta is a vector of censoring indicator
  # x is a nx2 matrix, with the 1st column being 1's and the 2nd column the score 
    
  ##--------------- define sub functions
  
  hazard.given.z=function(t, z, lamda, beta)
  {
    # t is survival time (either censored time or event time)
    # z is the latent variable indicating whether patient is in group 1 (high risk) or 0 (low risk)
    # lamda is the shape parameter for weibull distribution
    # beta=c(beta0, beta1), where beta0 is for the intercept and beta1 is for z
    
    return(lamda*t^(lamda-1)*exp(beta[1]+beta[2]*z))  
    
  }
  
  surv.given.z=function(t, z, lamda, beta)
  {
    # t is survival time (either censored time or event time)
    # z is the latent variable indicating whether patient is in group 1 or 0
    # lamda is the shape parameter for weibull distribution
    # beta=c(beta0, beta1)
    
    return(exp(-t^lamda*exp(beta[1]+beta[2]*z))) 
    
  }
  
  z.func=function(z, x, gamma)
  {
    # this gives probability of z, given x and parameter gamma
    # z can be 1, indicating high risk, or 0, indicating low risk
    # x=c(1, s), where s is the score value
    # gamma=c(gamma0, gamma1)
    
    return(exp(x%*%gamma)^z/(1+exp(x%*%gamma)))
    
  }
  
  
  like.given.z=function(y, delta, z, lamda, beta)
  {
    
    return(hazard.given.z(y, z, lamda, beta)^delta*surv.given.z(y, z, lamda, beta))
    
  }
  
  
  jointlike.func=function(y, delta, lamda, beta, z, x, gamma)
  {
    # this is the joint likelihood of y, delta, and z, where y has weibull and z has binomial distribution
    # z is the latent variable indicating whether patient is in group 1 or 0
    # y is survival time (either censored time or event time)
    # delta is censoring variable indicating whether patient with survival time y had an event (=1) or censored (=0)
    
    return(like.given.z(y, delta, z, lamda, beta)*z.func(z, x, gamma))
    
  }
  
  nQ.func=function(lam, b0, b1, g0, g1)
  {
    # this function is to calculate negative Q function for all the patients
    # the unknown variables are lam, b0, b1, g0, g1
    # fixed variables include y, delta, x and exp.p, where exp.p is expected probability 
    # of high risk given data and current parameter estimate
    
    tmp=cbind(1, exp.p)%*%c(b0,b1)
    
    tmp1=delta%*%(log(lam)+(lam-1)*log(y)+tmp)
    tmp2=-exp(b0)*y^lam%*%as.numeric(exp(b1)*exp.p+(1-exp.p))
    tmp3=as.numeric(exp.p)%*%as.numeric(log(z.func(1, x, c(g0, g1))))+as.numeric(1-exp.p)%*%as.numeric(log(1-z.func(1, x, c(g0, g1))))
    
    mysum=(tmp1+tmp2+tmp3)*(-1)
    return(mysum)
    
  }
  
  cat('\n....calculating parameter estimates using EM')
  
  ##extract initial values from the function arguments
  n=length(y)
  if (!is.null(initial.values))
  {
    z=initial.values[[1]]
    lamda=initial.values[[2]]
    beta=initial.values[[3]] 
    gamma=initial.values[[4]]   
  } else
  {
    z=rep(0,n)
    z[x[,2]>=quantile(x[,2], 0.75)]=1
    lamda=1
    beta=c(0, 0.5)
    gamma=c(0, 0.5)
        
  }  
  p.lamda=1
  p.beta=length(beta)
  p.gamma=length(gamma)
  p=p.lamda+p.beta+p.gamma #total dimension of unknown parameters
  
  k=0
  diff=.1
  loglike=sum(log(jointlike.func(y, delta, lamda, beta, 1, x, gamma)+jointlike.func(y, delta, lamda, beta, 0, x, gamma)))
  # both num and denom are vectors
  num=jointlike.func(y, delta, lamda, beta, 1, x, gamma)
  denom=jointlike.func(y, delta, lamda, beta, 1, x, gamma)+jointlike.func(y, delta, lamda, beta, 0, x, gamma)
  exp.p=num/denom 
  # a vector that tells the probability of high risk for each patient
  summary(exp.p)
  iter.results=c(k, lamda, beta, gamma, loglike, diff)
  iter.results
  em.error=NA
  
  while ((diff>error) & (k<max.iter))  
  {
    
    k=k+1; k
    #---update parameter estimates
    start=list(lam=as.numeric(lamda), 
               b0=as.numeric(beta[1]), 
               b1=as.numeric(beta[2]), 
               g0=as.numeric(gamma[1]), 
               g1=as.numeric(gamma[2]))
    res=try(mle(minuslogl=nQ.func, 
            start= start, 
            method = "BFGS"))
    if (class(res)=='try-error') {
      cat('\n......error has occured, so stopped EM!!!\n')
      em.error='yes'
      iter.results=rbind(iter.results, c(k, rep(9999,7)))
      k=max.iter
      
    } else {
      
      #summary(res)@coef
      vcov=res@vcov; vcov
      
      new.lamda=res@coef[1]
      new.beta=res@coef[2:3]
      new.gamma=res@coef[4:5]
      
      tmp=c(as.vector(abs(new.beta-beta)), as.vector(abs(new.gamma-gamma)), as.vector(abs(new.lamda-lamda)))
      tmp
      
      diff=max(tmp)
      diff
      
      #---assign new values to old
      lamda=new.lamda
      gamma=new.gamma
      beta=new.beta
      loglike=sum(log(jointlike.func(y, delta, lamda, beta, 1, x, gamma)+jointlike.func(y, delta, lamda, beta, 0, x, gamma)))
      
      #---update probability of high risk for every patient given current parameter estimates
      #---both num and denom are vectors
      num=jointlike.func(y, delta, lamda, beta, 1, x, gamma)
      denom=jointlike.func(y, delta, lamda, beta, 1, x, gamma)+jointlike.func(y, delta, lamda, beta, 0, x, gamma)
      exp.p=num/denom # a vector that tells the probability of high risk for each patient
      
      iter.results=rbind(iter.results, c(k, lamda, beta, gamma, loglike, diff))
            
    }  
    
  } #the end of the while loop

  #if after max.iter iterations, diff is still not <0.001 but it is <=0.01 then we say it converged, otherwise not 
  if ((k==max.iter) & (diff>0.01)) em.error='yes'
  
  ##above produced the final estimtes: lamda, gamma, beta, loglike, and 
  ##vcov as the final covariance matrix for the final estimates based on the nQ function;
  ##it also produced exp.p, which is from the last iteration
  
  colnames(iter.results)=c("iteration", "lamda", 
                           paste("b", 0:(p.beta-1), sep=''), 
                           paste("g", 0:(p.gamma-1), sep=''), "loglike", "diff")
  iter.results[1, 'diff']=NA
  
  if (is.na(em.error))
  {
    par.names=c('lamda', 'b0', 'b1', 'g0', 'g1')
    rownames(vcov)=colnames(vcov)=par.names
    
  } else if (em.error=='yes')
  {
    vcov=NA
  }
    
  out=list(em.iterations=iter.results, estimates.p=as.numeric(exp.p), nQ.vcov=vcov, em.error=em.error)
  return(out)
}


cov.func=function(em.results, y, delta, x)
{
  #################################
  # get covariance matrix using 
  # louis' formula with simulations.
  # need to use the EM estimates
  #################################
  
  cat('\n....calculating covariance matrix using Lous method and simulations')
  
  n=length(y); n
  estimates=(em.results$em.iterations)[nrow(em.results$em.iterations), c('lamda', 'b0', 'b1', 'g0', 'g1')]
  estimates
  exp.p=em.results$estimates.p
  vcov=em.results$nQ.vcov
  vcov
  
  lamda=estimates['lamda']; p.lamda=length(lamda)
  beta=estimates[c('b0', 'b1')]; p.beta=length(beta)
  gamma=estimates[c('g0', 'g1')]; p.gamma=length(gamma)
  b0=beta[1]
  b1=beta[2]
  g0=gamma[1]
  g1=gamma[2]
  
  ## now having got vcov as the covariance matrix of the nQ function 
  ## (minus complete loglikelihood), which is the first term in Lou's
  ## formula. We need to calculate the 2nd term now, that is, we need to 
  ## calculate the covariance matrix of the H function with simulations 
  ## (see Tanner and Wang, page 77 using simulation method)
  ## first generate z1,..., z10000, from bernoulli with probability of exp.p
  ## remember exp.p is a vector of n
  num.iter=10000
  z=matrix(NA, nrow=n, ncol=num.iter) 
  for (i in 1:n)
  {
    set.seed(i*1234)
    z[i,]=rbinom(n=num.iter, size=1, prob=exp.p[i])
  }
  
  ## Blow are partial derivatives of log complete likelihood.
  ## we want to simulate them wrt z, with probability exp.p from the last estimates
  p1.lamda=rep(NA, num.iter)
  p1.beta=matrix(NA, nrow=p.beta, ncol=num.iter)
  p1.gamma=matrix(NA, nrow=p.gamma, ncol=num.iter)
  p=p.lamda+p.beta+p.gamma
  
  pi.est=exp(x%*%gamma)/(1+exp(x%*%gamma))
  
  for (k in 1:num.iter)
  {
    m=cbind(rep(1, n), z[,k])
    p1.lamda[k]=delta%*%(1/lamda+log(y))-(y^lamda*log(y))%*%exp(m%*%beta)
    #first derivative of log of complete likelihood wrt lamda, 
    #notice difference between elementwise multiplication and vector operations
    
    p1.beta[,k]=matrix(delta-y^lamda*exp(m%*%beta), ncol=n)%*%m
    #first derivative of beta
    
    p1.gamma[,k]=matrix(z[,k]-pi.est, ncol=n)%*%x
    #first derivative of gamma
    
  }
  
  p1=rbind(p1.lamda, p1.beta, p1.gamma)  ## nrow=p, ncol=num.iter
  dim(p1)
  tmp=cov(t(p1)); tmp
  
  solve(vcov)
  solve(vcov)-tmp
  final.vcov=solve(solve(vcov)-tmp); final.vcov
  final.se=sqrt(diag(final.vcov)); final.se
  
  #estimate confidence interval using normal approximation 
  #(although lamda may not follow a normal approximation). 
  #Alternatively,can estimate log.lamda and then convert back to lamda
  #Then do hypothesis testing to test whether lamda=1, 
  #and all other parameters =0
  ll=estimates-1.96*final.se
  ul=estimates+1.96*final.se
  ll #lower limt
  ul
  test=(estimates-c(1,0,0,0,0))/final.se  
  #standardize and do the tests individually 
  test
  #P values using pnorm function, which uses lower tail as default
  pvalue=NULL
  pvalue[1]=2*(min(1-pnorm(test[1]), pnorm(test[1])))
  pvalue[2]=2*(min(1-pnorm(test[2]), pnorm(test[2])))
  pvalue[3]=2*(min(1-pnorm(test[3]), pnorm(test[3])))
  pvalue[4]=2*(min(1-pnorm(test[4]), pnorm(test[4])))
  pvalue[5]=2*(min(1-pnorm(test[5]), pnorm(test[5])))
  pvalue

  par.names=c('lamda', 'b0', 'b1', 'g0', 'g1')
  
  names(estimates)=par.names
  names(ll)=par.names
  names(ul)=par.names
  names(pvalue)=par.names
  names(final.se)=par.names
  rownames(final.vcov)=colnames(final.vcov)=par.names
  
  out=list(par.est=estimates, final.vcov=final.vcov, final.se=final.se,
           par.est.ll=ll, par.est.ul=ul, pvalue=pvalue, exp.p=exp.p)
  return(out)
}


greyzone.func=function(cov.results, y, delta, x, plot.logistic=T)
{
  ####################
  # find grey zone
  ####################
  
  cat('\n....finding grey zone in logistic curves\n')

  estimates=cov.results$par.est
  gamma=estimates[c('g0', 'g1')]
  final.se=cov.results$final.se
  final.vcov=cov.results$final.vcov
  exp.p=cov.results$exp.p
  
  #first estimated logistic regression curve and draw it as well as its CI
  pi.est=exp(x%*%gamma)/(1+exp(x%*%gamma)); summary(pi.est)
  score=x[,2]
  min.score=round(min(score),2)
  max.score=round(max(score),2)
  score.sim=seq(min.score,max.score,by=.01); length(score.sim)
  x.sim=cbind(rep(1,length(score.sim)), score.sim)
  #logistic curve based on score.sim (we want to get the curve at a more refined level than that based on score)
  pi.est.sim=exp(x.sim%*%gamma)/(1+exp(x.sim%*%gamma))
  
  #to find lower and upper bound for pi.est.sim curve, need to find estimated variance of g0+x*g1, or x.sim%*%gamma
  est.var=final.vcov[4,4]+score.sim^2*final.vcov[5,5]+2*score.sim*final.vcov[4,5]; est.var[1:10]
  est.sd=sqrt(est.var); est.sd[1:10]
  
  #95% CI of the curve of pi.est and using it to define a gray zone (third group)
  pi.est.ul.sim=exp(x.sim%*%gamma+1.96*est.sd)/(1+exp(x.sim%*%gamma+1.96*est.sd))
  pi.est.ll.sim=exp(x.sim%*%gamma-1.96*est.sd)/(1+exp(x.sim%*%gamma-1.96*est.sd))
  max.score.sim.lowrisk=score.sim[max((1:length(score.sim))[pi.est.ul.sim<0.5])]
  min.score.sim.highrisk=score.sim[min((1:length(score.sim))[pi.est.ll.sim>=0.5])]
  max.score.sim.lowrisk #
  min.score.sim.highrisk #
  p.low=exp(c(1,max.score.sim.lowrisk)%*%gamma)/(1+exp(c(1,max.score.sim.lowrisk)%*%gamma))
  p.high=exp(c(1,min.score.sim.highrisk)%*%gamma)/(1+exp(c(1,min.score.sim.highrisk)%*%gamma))
  as.numeric(p.low)  #
  as.numeric(p.high) #
  #so grey zone is >=max.score.sim.lowrisk and <min.score.sim.highrisk,
  #low risk is <max.score.sim.lowrisk, and high risk is >= min.score.sim.highrisk
  #however, greyzone doesn't have to exist, that is, max.score.sim.lowrisk and min.score.sim.highrisk can be NA
      
  if (plot.logistic)
  {
    plot(x=score, y=exp.p, cex=.6, xlim=c(-2, 3))
    points(x=score.sim, y=pi.est.sim, lty=1, col='blue', type='l')
    points(x=score.sim, y=pi.est.ul.sim, lty=2, col='blue', type='l')
    points(x=score.sim, y=pi.est.ll.sim, lty=2, col='blue', type='l')
    abline(h=0.5)
    abline(v=c(max.score.sim.lowrisk, min.score.sim.highrisk), col='green', lty=2)
    
  }
  
  out=list(greyzone.ll=max.score.sim.lowrisk, greyzone.ul=min.score.sim.highrisk,
           score=score, score.sim=score.sim, pi.est.sim=pi.est.sim, exp.p=exp.p,
           pi.est.ul.sim=pi.est.ul.sim, pi.est.ll.sim=pi.est.ll.sim) 
    
}


cox.summary=function(stime, sind, var)
{
  
  p=length(unique(var))-1
  
  #calculate r2  
  data=data.frame(stime, sind, var)
  n=length(stime)
  model0 <- coxph(Surv(stime, sind)~1)
  model1 <- coxph(Surv(stime, sind)~as.factor(var))
  f0 <- rep(0,n)
  f1 <- predict(model1, newdata=data)
  Surv.res <- Surv(stime, sind)
  r2.oxs=survAUC::OXS(Surv.res, f1, f0)
  r2.nagelk=survAUC::Nagelk(Surv.res, f1, f0)
  r2.xo=survAUC::XO(Surv.res, f1, f0)
  
  res=summary(model1); res
  pvalue=as.numeric(res$logtest[3])    
  aic=-2*res$loglik[2]+2*p
  
  #considering ties in x
  c1=Hmisc::rcorrcens(Surv(stime, sind) ~ var); c1
  c1=1-c1[1,1];c1
  #ignores ties in x
  c3=Hmisc::rcorr.cens(var, Surv(stime, sind), outx=T)
  c3=1-c3['C Index'];c3
  
  return(out=c(p=pvalue, AIC=aic, r2.oxs=r2.oxs, r2.xo=r2.xo, r2.nagelk=r2.nagelk, 
               c.ties=c1))    
  
}


bestcut2=function(data=data, stime="surtim", sind="surind", var="score", leave=20)
{
  # this function finds the best cutoff value for a continuous score given survival data to divide patients into 2 groups
  #leave is the minimum number of patients in each patient group
  
  cat('\nLooking for best cutoff to divide subjects into 2 groups...')
  n=nrow(data)
  ord=order(data[,var]); ord
  score=data[ord,var]; score
  
  pvalue=chisq.stat=rep(NA, n)
  for (i in (leave+1):(n-leave+1))
  {  
    ## i is the last position of the 1st cluster
    
    cluster=rep(0,n)
    cluster[data[,var]>=score[i]]=1
    res=survdiff(Surv(data[,stime], data[,sind])~cluster)
    pvalue[i]=1-pchisq(res$chisq, df=1)
    chisq.stat[i]=res$chisq
    
  }
  summary(-log10(pvalue))
    
  ##### which cutoff gives the smallest p-value: best.i gives the answer
  best.i=which.max(chisq.stat)
  best.i
  cutoff=score[best.i]
  cat('\nbest cutoff =', cutoff)
  
  cluster = rep(0,n)
  cluster[data[,var]>=cutoff]=1
  
  survdata=data.frame(data, pvalue, chisq.stat, bestcut2=cluster, cutoff=rep(cutoff,n))

### are there pvalues significant after Bonferroni?

  totsig=sum(!is.na(pvalue))
  siglevel=-log10(.05/totsig)  # the significance level after Bonferroni correction
  siglevel
  maxp=max(as.vector(-log10(pvalue)),na.rm=T)
  maxp
  cat('\nThere are p values significant after Bonferroni correction:', maxp>siglevel)
  cat("\nP-value at best cutoff: ", pvalue[best.i], "\n\n")
  
  return(survdata)
  
}

