# bivariate standard normal log-likelihood for one observation
# input:
# low the vector of lower limits of length n.
# upp the vector of upper limits of length n.
# r the correlation parameter
# output:
# bivariate standard normal log-likelihood for one observation
bivlik<-function(low,upp,r)
{ rmat<-matrix(c(1,r,r,1),2,2)
  prob<-pmvnorm(lower=low,upper=upp,mean=rep(0,2),corr=rmat)[1]
  log(prob)
}


# the mean values of the univariate marginal distribution 
# corresonding to the used link function 
# input:
# x the matix of the covariates 
# b the vector with the regression coefficients
# link has three options: 1. "log", 2. "logit". 3. "probit"
# output:
# the mean values of the univariate marginal distribution 
linked.mu<-function(x,b,link)
{ if(link=="log")
    { mu<-exp(x %*% b)
    }
    else
    { if(link=="logit")
      { expnu<-exp(x %*% b)
        mu<-expnu/(1+expnu)
      }
    else
    { # link=probit
      mu<-pnorm(x %*% b) }
    }
   mu
}

# Density  of the univariate marginal distribution
# input:
# y the vector of (non-negative integer) quantiles.
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# invgam the inverse of parameter  gamma of negative binomial distribution.
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
# output:
# the density of the univariate marginal  distribution
dmargmodel<-function(y,mu,gam,invgam,margmodel)
{ if(margmodel=="poisson")
    { dpois(y,mu)
    }
    else
    { if(margmodel=="bernoulli")
      { dbinom(y,size=1,prob=mu) }
    else
    { if(margmodel=="nb1")
      { dnbinom(y,prob=1/(1+gam),size=mu*invgam) }
    else
    { # margmodel=="nb2"
      dnbinom(y,size=invgam,mu=mu) }}}
}

# CDF  of the univariate marginal distribution
# input:
# y the vector of (non-negative integer) quantiles.
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# invgam the inverse of parameter  gamma of negative binomial distribution.
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
# output:
# the cdf of the univariate marginal  distribution
pmargmodel<-function(y,mu,gam,invgam,margmodel)
{ if(margmodel=="poisson")
    { ppois(y,mu)
    }
    else
    { if(margmodel=="bernoulli")
      { pbinom(y,size=1,prob=mu) }
    else
    { if(margmodel=="nb1")
      { pnbinom(y,prob=1/(1+gam),size=mu*invgam) }
    else
    { # margmodel=="nb2"
      pnbinom(y,size=invgam,mu=mu) }}}
}

# quantile  of the univariate marginal distribution
# input:
# y the vector of probabilities
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# invgam the inverse of parameter  gamma of negative binomial distribution.
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
# output:
# the quantile of the univariate marginal  distribution
qmargmodel<-function(y,mu,gam,invgam,margmodel)
{ if(margmodel=="poisson")
    { qpois(y,mu)
    }
    else
    { if(margmodel=="bernoulli")
      { qbinom(y,size=1,prob=mu) }
    else
    { if(margmodel=="nb1")
      { qnbinom(y,prob=1/(1+gam),size=mu*invgam) }
    else
    { # margmodel=="nb2"
      qnbinom(y,size=invgam,mu=mu) }}}
}

# negative univariate logikelihood assuming independence within clusters
# input:
# param the vector of regression and not regression parameters
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# negative univariate logikelihood assuming independence within clusters
marglik<-function(param,xdat,ydat,margmodel,link)
{ p<-dim(xdat)[2]
  b<-param[1:p]
  if(margmodel=="nb1" | margmodel=="nb2")
  { gam<-param[p+1]
    invgam<-1/gam
  }
  #else  gam<-invgam<-0
  mu<-linked.mu(as.matrix(xdat),b,link)
  -sum(log(dmargmodel(ydat,mu,gam,invgam,margmodel)))
}
 
 
# Independent estimating equations for binary, Poisson or 
# negative binomial regression.
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output: A list containing the following components:
# coef the vector with the ML estimated regression parameters
# gam the ML estimate of gamma parameter
iee<-function(xdat,ydat,margmodel,link="log")
{ #if(margmodel=="bernoulli")  family=binomial else family=poisson
  if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  if(margmodel=="poisson")
  { uni<-glm(ydat ~ xdat[,-1],family =poisson(link="log"))
    res<-as.vector(uni$coef)
    list(reg=res)
  } else {
  if(margmodel=="bernoulli")
  { if(link=="probit") 
    { uni<-glm(ydat ~ xdat[,-1],family =binomial(link="probit"))
    } else {
    uni<-glm(ydat ~ xdat[,-1],family =binomial(link="logit")) }     
    res<-as.vector(uni$coef)
    list(reg=res)
  } else
  { p<-dim(xdat)[2]
    uni<-nlm(marglik,c(rep(0,p),1),margmodel=margmodel,
    link=link,xdat=xdat,ydat=ydat,iterlim=1000)
    res1<-uni$e[1:p]
    res2<-uni$e[p+1]
    list(reg=res1,gam=res2) }}
}

# corralation matrix
# input:
# d the dimension
# r a vector with correaltion parameters
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# output:
# the correlation matrix
cormat<-function(d,r,pairs,corstr)
{ rmat<-matrix(1,d,d)
  lr<-nrow(pairs)
  if(length(r)==1) r<-rep(r,d^2)
  if(corstr=="exch" | corstr=="unstr")
  {
  for(j in 1:lr)
  { rmat[pairs[j,1],pairs[j,2]]<-r[j]
    rmat[pairs[j,2],pairs[j,1]]<-r[j]
  }
  } else {
  for(j in 1:lr)
  { temp<-pairs[j,2]-pairs[j,1]
    rmat[pairs[j,1],pairs[j,2]]<-r[j]^temp
    rmat[pairs[j,2],pairs[j,1]]<-r[j]^temp
  }
  }  
  rmat
}

# Bivariate composite likelihood for multivariate normal copula with Poisson, 
# binary, or negative binomial regression.
# input:
# r the vector of normal copula parameters
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# output:
# negative bivariate composite likelihood for multivariate normal copula 
# with Poisson, binary, or negative binomial regression.
bcl<-function(r,b,gam,xdat,ydat,id,tvec,margmodel,corstr,link="log")
{ s<-0
  if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit" 
  if(margmodel=="nb1" | margmodel=="nb2") invgam<-1/gam else invgam<-NULL
  uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  tvec<-id.time(tvec,d)
  n<-1:dim(xdat)[1]
  pairs<-maxpairs(d)
  rmat<-cormat(maxd,r,pairs,corstr)            
  c1<-(r[1] > 1 | r[1] < (-1/(maxd-1))) 
  c2<-(r[1] > 1 | r[1] < -1)            
  c3<-(det(rmat)<0 | min(rmat)< -1 | max(rmat)>1)
  if((corstr=="exch" & c1) | (corstr=="ar" & c2 ) | (corstr=="unstr" & c3)) {s<-1e10}
  else {  
  for(i in uid)
  { cases<-id==i
    irow=n[cases]
    yi<-ydat[irow]
    ti<-tvec[irow]
    newyi<-rep(NA,maxd)
    newmui<-rep(NA,maxd)
    newyi[ti]<-yi
    x<-xdat[cases,]
    mui<-linked.mu(x,b,link)
    newmui[ti]<-mui
    vlow<-pmargmodel(newyi-1,newmui,gam,invgam,margmodel)
    tem<-dmargmodel(newyi,newmui,gam,invgam,margmodel)
    vupp<-vlow+tem
    zlow=qnorm(vlow)
    zupp=qnorm(vupp)
    for(j in 1:dim(pairs)[1])
    { k1<-pairs[j,][1]
      k2<-pairs[j,][2]
      if((sum(k1==ti)==1) & (sum(k2==ti)==1))
      { if(corstr=="exch")
        { s<-s +bivlik(zlow[c(k1,k2)],zupp[c(k1,k2)],r) }
        else
        { if(corstr=="ar")
        { s<-s +bivlik(zlow[c(k1,k2)],zupp[c(k1,k2)],r^(k2-k1)) }
        else
        { # corstr=="unstr"
          s<-s +bivlik(zlow[c(k1,k2)],zupp[c(k1,k2)],r[j])
        }}}
      else { s<-s }
     }
  }
 -s
 }
}

# optimization routine for composite likelihood for MVN copula
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# output: A list containing the following components:
# minimum the value of the estimated minimum of CL1 for MVN copula
# estimate the CL1 estimates
# gradient the gradient at the estimated minimum of CL1
# code an integer indicating why the optimization process terminated, see nlm.
cl1<-function(b,gam,xdat,ydat,id,tvec,margmodel,corstr,link="log")
{ if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  d<-id.size(id)
  pairs<-maxpairs(d)
  if(corstr=="unstr")
  { nom<-dim(pairs)[1]
    nlm(bcl,rep(0.1,nom),b,gam,xdat,ydat,id,tvec,margmodel,corstr,link)
  }
  else
  {nlm(bcl,0.1,b,gam,xdat,ydat,id,tvec,margmodel,corstr,link)}
}



# the dimension of each id
# input: the vector with the id
# id the vector with the id
# output:
# a vector with the dimension of each id
id.size<-function(id)
{ d<-NULL
  uid<-unique(id)
  for(i in uid)
  { di<-sum(id==i)
    d<-c(d,di) }
 d
}


# the transformed time points, i.e, 1,2,3... for each subject
# input:
# tvec: untrasformed time points
# d the dimension for all subjects
# output:
# the transformed time points, i.e, 1,2,3... for each subject
id.time<-function(tvec,d)
{ maxd<-max(d)
  ut<-sort(unique(tvec))
  newtvec<-rep(NA,length(tvec))
  for(j in 1:maxd)
  { newtvec[tvec==ut[j]]<-j }
  newtvec
}

# the maximum number of bivariate pairs
# input:
# d the dimension for all subjects
# output:
# the maximum number of bivariate pairs
maxpairs<-function(d)
{ pairs<-NULL
  maxd<-max(d)
  for(id1 in 1:(maxd-1))
  { for(id2 in (id1+1):maxd)
  { pairs<-rbind(pairs,c(id1,id2)) } }
  pairs
}


# derivative of the marginal loglikelihood with respect to nu
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# ub the truncation value
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson,
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function,
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# the vector with the derivatives of the margmodel loglikelihood with respect to nu
derlik.nu<-function(mu,gam,invgam,ub,margmodel,link)
{ if(link=="probit")
  { if(mu==1 & is.finite(mu)){mu<-0.9999}
    nu<-qnorm(mu)
    (0:1-mu)/mu/(1-mu)*dnorm(nu)
  }
  else
  { if(margmodel=="nb1")
    { j<-0:(ub-1)
      s<-c(0,cumsum(1/(mu+gam*j)))
      (s-invgam*log(1+gam))*mu
    }
  else
  { if(margmodel=="nb2")
    { pr<-1/(mu*gam+1)
      (0:ub-mu)*pr
    }
  else { 0:ub-mu }}}
}





# derivative of the marginal loglikelihood with respect to nu
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# y  the value of a non-negative integer quantile
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson,
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function,
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# the derivative of the margmodel loglikelihood with respect to nu
iderlik.nu<-function(mu,gam,invgam,y,margmodel,link)
{ if(link=="probit")
  { if(mu==1 & is.finite(mu)){mu<-0.9999}
    nu<-qnorm(mu)
    (y-mu)/mu/(1-mu)*dnorm(nu)
  }
  else
  { if(margmodel=="nb1")
    { s<-0
      if(y>0)
      { j<-0:(y-1)
        s<-sum(1/(mu+gam*j))
      }
      (s-invgam*log(1+gam))*mu
    }
  else
  { if(margmodel=="nb2")
    { pr<-1/(mu*gam+1)
      (y-mu)*pr
    }
  else { y-mu }}}
}

# derivative of the NB loglikelihood with respect to gamma
# input:
# gam the gamma parameter
# invgam the inverse of gamma parameter
# mu the mean parameter
# ub the truncation value
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the vector with the derivatives of the NB loglikelihood with respect to gamma
derlik.gam<-function(mu,gam,invgam,ub,margmodel)
{ j<-0:(ub-1)
  if(margmodel=="nb1")
  { s<-c(0,cumsum(j/(mu+gam*j)))
    s+invgam*invgam*mu*log(1+gam)-(0:ub+invgam*mu)/(1+gam)
  }
  else
  { #if(margmodel=="nb2")
    pr<-1/(mu*gam+1)
    s<-c(0,cumsum(j/(1+j*gam)))
    s-log(pr)/(gam*gam)-(0:ub+invgam)*mu*pr
  }
}

# derivative of the NB loglikelihood with respect to gamma
# input:
# gam the gamma parameter
# invgam the inverse of gamma parameter
# mu the mean parameter
# y  the value of a non-negative integer quantile
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the derivative of the NB loglikelihood with respect to gamma
iderlik.gam<-function(mu,gam,invgam,y,margmodel)
{ s<-0
  if(margmodel=="nb1")
  { if(y>0)
  { j<-0:(y-1)
    s<-sum(j/(mu+gam*j))
  }
  s+invgam*invgam*mu*log(1+gam)-(y+invgam*mu)/(1+gam)
  }
  else
  { #if(margmodel=="nb2")
    if(y>0)
  { j<-0:(y-1)
    s<-sum(j/(1+gam*j))
  }
  pr<-1/(mu*gam+1)
  s-log(pr)/(gam*gam)-(y+invgam)*mu*pr
  }
}



# minus expectation of the second derivative of the marginal loglikelihood
# with resect to nu
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# u the univariate cdfs
# ub the truncation value
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the vector with the minus expectations of the margmodel loglikelihood
# with respect to nu
fisher.nu<-function(mu,gam,invgam,u,ub,margmodel,link)
{ if(link=="log" & margmodel=="poisson")
    { mu }
    else
    { if(link=="logit")
      { mu*(1-mu) }
    else
    { if(margmodel=="nb1")
      { j<-0:ub
        s1<-sum(1/(mu+j*gam)/(mu+j*gam)*(1-u))
        s2<-sum(1/(mu+j*gam)*(1-u))
        (mu*s1-s2+invgam*log(1+gam))*mu
      }
    else
    { if(margmodel=="nb2")
      { pr<-1/(mu*gam+1)
        mu*pr
      }
    else
    { # link=="probit"
      if(mu==1 & is.finite(mu)){mu<-0.9999}
      nu<-qnorm(mu)
      1/mu/(1-mu)*dnorm(nu)^2}
    }}}
}




# minus expectation of the second derivative of the marginal NB loglikelihood
# with resect to gamma
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# u the univariate cdfs
# ub the truncation value
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the vector with the minus expectations of the margmodel loglikelihood
# with respect to gamma
fisher.gam<-function(mu,gam,invgam,u,ub,margmodel)
{ j<-0:ub
  if(margmodel=="nb1")
  { pr<-1/(gam+1)
    s<-sum((j/(mu+j*gam))^2*(1-u))
    s+2*invgam*invgam*invgam*log(1+gam)*mu-2*invgam*invgam*mu*pr-mu/gam*pr
  }
  else
  { #if(margmodel=="nb2")
    s<-sum((invgam+j)^(-2)*(1-u))
    invgam^4*(s-gam*mu/(mu+invgam))
  }
}






# minus expectation of the second derivative of the marginal loglikelihood
# with resect to nu and gamma
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# u the univariate cdfs
# ub the truncation value
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the vector with the minus expectations of the NB loglikelihood
# with respect to nu and gamma
fisher.nu.gam<-function(mu,gam,invgam,u,ub,margmodel)
{ if(margmodel=="nb1")
  { pr<-1/(gam+1)
    j<-0:ub
    s<-sum(j/(mu+j*gam)/(mu+j*gam)*(1-u))
    (s-invgam*invgam*log(1+gam)+invgam*pr)*mu
  }
  else {0}
}

# Calculating the truncation value for the univariate distribution
# input:
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
# output:
# the truncation value--upper bound
truncation<-function(mu,gam,margmodel)
{ if(margmodel=="poisson")
    { ub<-round(max(10,mu+7*sqrt(mu),na.rm=T))
    }
    else
    { if(margmodel=="bernoulli") ub<-1
    else
    { if(margmodel=="nb1")
      { pr<-1/(gam+1)
        v<-mu/pr
        ub<-round(max(10,mu+10*sqrt(v),na.rm=T))
       }
    else
    { pr<-1/(mu*gam+1)
      v<-mu/pr
      ub<-round(max(10,mu+7*sqrt(v),na.rm=T))
    }}}
  ub
}



# approximation of bivariate normal cdf (Johnson&Kotz, 1972)
# For rho<=0.4 the series truncated   at rho^3.
# For larger rho truncation at rho^5.
# input:
# r the normal copula parameter
approxbvncdf<-function(r,x1,x2,x1s,x2s,x1c,x2c,x1f,x2f,t1,t2)
{ r2=r*r; r3=r*r2; r4=r2*r2; r5=r4*r
  tem3=r+r2*outer(x1,x2)/2+r3*outer(x1s-1,x2s-1)/6
  if(r>0.4)
  { tem4=tem3+r4*outer(x1c-3*x1, x2c-3*x2)/24
    tem5=tem4+r5*outer(x1f-6*x1s+3, x2f-6*x2s+3)/120
    pr.a5=t1+t2*tem5
    pr.a5
  } else { pr.a3=t1+t2*tem3; pr.a3}
}


# covariance matrix of the scores Omega_i
# input:
# scnu the matrix of the score functions with respect to nu
# scgam the matrix of the score functions with respect to gam
# index the bivariate pair
# pmf the matrix of rectangle probabilities
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
scoreCov<-function(scnu,scgam,pmf,index,margmodel)
{ j1<-index[1]
  j2<-index[2]
  cov11<-t(scnu[,j1])%*%pmf%*%scnu[,j2]
  if(margmodel=="bernoulli" | margmodel=="poisson")
  { cov11 }
  else
  { cov12<-t(scnu[,j1])%*%pmf%*%scgam[,j2]
    cov21<-t(scgam[,j1])%*%pmf%*%scnu[,j2]
    cov22<-t(scgam[,j1])%*%pmf%*%scgam[,j2]
    matrix(c(cov11,cov12,cov21,cov22),2,2)
  }
}



# select the present column and lines in Omega, Delta and X matrices
# for unbalanced data
# input:
# tvec a vector of the time for an individual
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
subselect<-function(tvec,margmodel)
{ if(margmodel=="bernoulli" | margmodel=="poisson") sel<-tvec
  else
  { sel<-NULL
    for(i in 1:length(tvec))
    { tm<-2*tvec[i]
      k<-c(tm-1,tm)
      sel<-c(sel,k)
    }}
  sel
}



# weight matrix fixed at values from the CL1 estimator
# input:
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# rh the vector with CL1 estimates
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# output: A list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
weightMat<-function(b,gam,rh,xdat,ydat,id,tvec,margmodel,corstr,link="log")
{ if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  if(margmodel=="nb1" | margmodel=="nb2") invgam<-1/gam else invgam<-NULL
  uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  q<-length(gam)
  lid<-length(uid)
  dim<-dim(xdat)
  n<-1:dim[1]
  p<-dim[2]
  pairs<-maxpairs(d)
  tvec<-id.time(tvec,d)
  omega<-array(NA,c(maxd*(1+q),maxd*(1+q),lid))
  X<-array(NA,c(p+q,maxd*(1+q),lid))
  delta<-array(NA,c(maxd*(1+q),maxd*(1+q),lid))
  dom<-maxd*(1+q)
  pos<-seq(1,dom-1,by=2)    #not used for binary
  m<-0
  for(i in uid)
  { m<-m+1
    cases<-id==i
    irow=n[cases]
    newx<-matrix(NA,maxd,p)
    newyi<-rep(NA,maxd)
    newmui<-rep(NA,maxd)
    ti<-tvec[irow]
    x<-xdat[cases,]
    mui<-linked.mu(x,b,link)
    newx[ti,]<-x
    newmui[ti]<-mui
    ub<-truncation(newmui,gam,margmodel)
    du<-scnu<-scgam<-matrix(NA,1+ub,maxd)
    for(j in ti)
    { du[,j]<-dmargmodel(0:ub,newmui[j],gam,invgam,margmodel)
      scnu[,j]<-derlik.nu(newmui[j],gam,invgam,ub,margmodel,link)
      if(margmodel=="nb1" | margmodel=="nb2")
      { scgam[,j]<-derlik.gam(newmui[j],gam,invgam,ub,margmodel) }
    }
    u<-apply(du,2,cumsum)
    z<-qnorm(u)
    z[is.nan(z)]<-7
    z[z>4&margmodel=="bernoulli"]<-7
    z<-rbind(-7,z)
    pz<-pnorm(z)
    dz<-dnorm(z)
    zs<-z*z
    zc<-z*zs
    zf<-zs*zs
    xi<-NULL
    if(margmodel=="bernoulli" | margmodel=="poisson")
    { diagonali<-rep(NA,maxd)
    } else {
    diagonali<-array(NA,c(2,2,maxd)) }
    for(j in 1:maxd)
    { f1<-fisher.nu(newmui[j],gam,invgam,u[,j],ub,margmodel,link)
      if(margmodel=="nb1" | margmodel=="nb2")
      { temp<-cbind(newx[j,],0)
        xi<-cbind(xi,rbind(temp,c(0,1)))
        f2<-fisher.gam(newmui[j],gam,invgam,u[,j],ub,margmodel)
        f3<-fisher.nu.gam(newmui[j],gam,invgam,u[,j],ub,margmodel)
        diagonali[,,j]<-matrix(c(f1,f3,f3,f2),2,2)
      }
      else
      { temp<-newx[j,]
        xi<-cbind(xi,temp)
        diagonali[j]<-f1
      }
    }
    if(margmodel=="bernoulli" | margmodel=="poisson")
    { deltai<-diag(diagonali)
      offi<-rep(NA,dim(pairs)[1])
    } else {
      deltai<-matrix(0,dom,dom)
      minus<-0
      for(j in pos)
      { deltai[j:(j+1),j:(j+1)]<-diagonali[,,j-minus]
        minus<-minus+1
      }
      offi<-array(NA,c(2,2,dim(pairs)[1]))
    }
    for(k in 1:dim(pairs)[1])
    { k1<-pairs[k,][1]
      k2<-pairs[k,][2]
      if((sum(k1==ti)==1) & (sum(k2==ti)==1))
      { x1=z[,k1]; x2=z[,k2]
        x1s=zs[,k1]; x2s=zs[,k2]
        x1c=zc[,k1]; x1f=zf[,k1]
        x2c=zc[,k2]; x2f=zf[,k2]
        t1=outer(pz[,k1],pz[,k2])
        t2=outer(dz[,k1],dz[,k2])
        if(corstr=="exch")
        { cdf<-approxbvncdf(rh,x1,x2,x1s,x2s,x1c,x2c,x1f,x2f,t1,t2) }
        else
        { if(corstr=="ar")
        { cdf<-approxbvncdf(rh^(k2-k1),x1,x2,x1s,x2s,x1c,x2c,x1f,x2f,t1,t2) }
        else
        { # corstr=="unstr"
          cdf<-approxbvncdf(rh[k],x1,x2,x1s,x2s,x1c,x2c,x1f,x2f,t1,t2)
        }}
        cdf1=apply(cdf,2,diff)
        pmf=apply(t(cdf1),2,diff)
        pmf=t(pmf)
        if(margmodel=="bernoulli" | margmodel=="poisson")
        {offi[k]<-scoreCov(scnu,scgam,pmf,pairs[k,],margmodel)}
        else {offi[,,k]<-scoreCov(scnu,scgam,pmf,pairs[k,],margmodel)}
      }
    }
    omegai<-deltai
    if(margmodel=="bernoulli" | margmodel=="poisson")
    { for(j in 1:dim(pairs)[1])
      { omegai[pairs[j,1],pairs[j,2]]<-offi[j]
        omegai[pairs[j,2],pairs[j,1]]<-offi[j]
      }}
    else
    { ch1<-0
      ch2<-0
      for(j in 1:(maxd-1))
      { for(r in pos[-(1:j)])
        { omegai[(j+ch1):(j+1+ch1),r:(r+1)]<-offi[,,(j+ch2-ch1)]
          omegai[r:(r+1),(j+ch1):(j+1+ch1)]<-t(offi[,,(j+ch2-ch1)])
          ch2<-ch2+1
        }
      ch1<-ch1+1
    }}
    X[,,m]<-xi
    delta[,,m]<-deltai
    omega[,,m]<-omegai
    }
    list(omega=omega,X=X,delta=delta)
}


# the weigted scores equations
# input:
# param the vector of regression and not regression parameters
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# WtScMat is a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# output
# the weigted scores equations
wtsc<-function(param,WtScMat,xdat,ydat,id,tvec,margmodel,link="log")
{ if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  dim<-dim(xdat)
  n<-1:dim[1]
  p<-dim[2]
  tvec<-id.time(tvec,d)
  b<-param[1:p]
  if(p<length(param)) {gam<-param[p+1]; invgam<-1/gam }
  g<-0
  m<-0
  for(i in uid)
  { cases<-id==i
    irow=n[cases]
    m<-m+1
    x<-xdat[cases,]
    mu<-linked.mu(x,b,link)
    ub<-truncation(mu,gam,margmodel)
    scnu<-scgam<-matrix(NA,ub+1,d[m])
    for(j in 1:d[m])
    { scnu[,j]<-derlik.nu(mu[j],gam,invgam,ub,margmodel,link)
      if(margmodel=="nb1" | margmodel=="nb2")
      { scgam[,j]<-derlik.gam(mu[j],gam,invgam,ub,margmodel) }
    }
    y<-ydat[irow]
    sci<-NULL
    for(j in 1:d[m])
    { if(y[j]>ub)
      { scnui<-iderlik.nu(mu[j],gam,invgam,y[j],margmodel,link)
        if(margmodel=="nb1" | margmodel=="nb2")
        { scgami<-iderlik.gam(mu[j],gam,invgam,y[j],margmodel)
        } else {
        scgami<-NULL}
      }
      else {
      scnui<-scnu[y[j]+1,j]
      if(margmodel=="nb1" | margmodel=="nb2")
      { scgami<-scgam[y[j]+1,j]
      } else {
        scgami<-NULL}}
        sci<-c(sci,c(scnui,scgami))
    }
    ti<-tvec[irow]
    seli<-subselect(ti,margmodel)
    Xi<-WtScMat$X[,,m]
    Xi<-Xi[,seli]
    deltai<-WtScMat$delta[,,m]
    deltai<-deltai[seli,seli]
    omegai<-WtScMat$omega[,,m]
    omegai<-omegai[seli,seli]
    gi<-Xi%*%t(deltai)%*%solve(omegai,sci)
    g<-g+gi
    }
    g
}


# solving the weigted scores equations
# input:
# start the starting values (IEE estimates) for the vector of
# regression and not regression parameters
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# WtScMat is a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# output:
# the weighted scores estimates
solvewtsc<-function(start,WtScMat,xdat,ydat,id,tvec,margmodel,link="log")
{ multiroot(f=wtsc,start,atol=1e-4,rtol=1e-4,ctol=1e-4,
WtScMat=WtScMat,xdat=xdat,ydat=ydat,id=id,tvec=tvec,margmodel=margmodel,link=link) }



# inverse Godambe matrix with Delta and Omega evaluated at IEE estimator
# input:
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# rh the vector with CL1 estimates
# WtScMat a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# the inverse Godambe matrix
godambe<-function(param,WtScMat,xdat,ydat,id,tvec,margmodel,link="log")
{ if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit" 
  uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  dim<-dim(xdat)
  n<-1:dim[1]
  p<-dim[2]
  tvec<-id.time(tvec,d)
  b<-param[1:p]
  if(margmodel=="nb1" | margmodel=="nb2")
  { gam<-param[p+1]
    invgam<-1/gam
  } else {
  gam<-invgam<-NULL}
  v<-v1<-v2<-v3<-0
  fv<-fv1<-fv2<-fv3<-0
  m<-0
  for(i in uid)
  { m<-m+1
    cases<-id==i
    irow=n[cases]
    newx<-matrix(NA,max(d),length(b))
    newy<-rep(NA,maxd)
    newmui<-rep(NA,maxd)
    ti<-tvec[irow]
    x<-xdat[cases,]
    mui<-linked.mu(x,b,link)
    newx[ti,]<-x
    newmui[ti]<-mui
    ub<-truncation(newmui,gam,margmodel)
    scnu<-scgam<-matrix(NA,1+ub,maxd)
    for(j in ti)
    { scnu[,j]<-derlik.nu(newmui[j],gam,invgam,ub,margmodel,link)
      if(margmodel=="nb1" | margmodel=="nb2")
      { scgam[,j]<-derlik.gam(newmui[j],gam,invgam,ub,margmodel) }
    }
    y<-ydat[irow]
    newy[ti]<-y
    sci<-NULL
    for(j in 1:maxd)
    { if(sum(newy[j]>ub,na.rm=T)==1)
      { scnui<-iderlik.nu(newmui[j],gam,invgam,newy[j],margmodel,link)
        if(margmodel=="nb1" | margmodel=="nb2")
        { scgami<-iderlik.gam(newmui[j],gam,invgam,newy[j],margmodel)
        } else {
        scgami<-NULL}
      }
      else {
      scnui<-scnu[newy[j]+1,j]
      if(margmodel=="nb1" | margmodel=="nb2")
      { scgami<-scgam[newy[j]+1,j]
      } else {
        scgami<-NULL}}
        sci<-c(sci,c(scnui,scgami))
    }
    sci<-sci[!is.na(sci)]
    seli<-subselect(ti,margmodel)
    Xi<-WtScMat$X[,,m]
    Xi<-Xi[,seli]
    fdeltai<-WtScMat$delta[,,m]
    fdeltai<-fdeltai[seli,seli]
    fomegai<-WtScMat$omega[,,m]
    fomegai<-fomegai[seli,seli]
    fci<-fdeltai%*%t(Xi)
    nrhs1=ncol(fdeltai)
    nrhs2=ncol(fci)
    tem.lin=solve(fomegai,cbind(fdeltai,fci))
    finvai=t(tem.lin[,1:nrhs1])
    xi.invai=Xi %*% finvai
    fai<-solve(finvai)
    tem=solve(t(fai),t(Xi))
    fv1i<-xi.invai %*% fci
    fv2i<-xi.invai %*% (sci%*% t(sci))  %*% tem
    fv3i<-t(fci) %*% tem
    fv1<-fv1+fv1i
    fv2<-fv2+fv2i
    fv3<-fv3+fv3i
    }
    solve(fv1,fv2)%*%solve(fv3)
}



# the weighted scores wrapper function: handles all the steps in the weighted scores
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# rh the vector with CL1 estimates
# WtScMat a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# iprint indicates printing of some intermediate results, default FALSE
wtsc.wrapper<-function(xdat,ydat,id,tvec,margmodel,corstr,link="log",iprint=FALSE)
{ if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  i.est<-iee(xdat,ydat,margmodel,link)
  if(iprint)
  { cat("\niest: IEE estimates\n")
    print(c(i.est$reg,i.est$gam))
  }
  est.rho<-cl1(b=i.est$reg,gam=i.est$gam,xdat,ydat,id,tvec,margmodel,corstr,link)
  if(iprint)
  { cat("\nest.rho: CL1 estimates\n")
    print(est.rho$e)
    cat("\nest.rho: CL1 likelihood\n")
    print(-est.rho$m)
  }
  WtScMat<-weightMat(b=i.est$reg,gam=i.est$gam,rh=est.rho$e,
           xdat,ydat,id,tvec,margmodel,corstr,link)
  ws<-solvewtsc(start=c(i.est$reg,i.est$gam),WtScMat,xdat,ydat,id,
  tvec,margmodel,link)
  if(iprint)
  { cat("ws=parameter estimates\n")
    print(ws$r)
  }
  acov<-godambe(ws$r,WtScMat,xdat,ydat,id,tvec,margmodel,link)
  se<-sqrt(diag(acov))
  if(iprint)
  { cat("\nacov: inverse Godambe matrix with W based on first-stage wt matrices\n")
    print(acov)
    cat("\nse: robust standard errors\n")
    print(se)
    res<-round(cbind(ws$r,se),3)
    cat("\nres: Weighted scores estimates and standard errors\n")
    print(res)
  }
  list(IEEest=c(i.est$reg,i.est$gam),CL1est=est.rho$e,CL1lik=-est.rho$m,WSest=ws$r, asympcov=acov)
}


