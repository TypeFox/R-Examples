# v0.0.1 created on Feb. 10, 2009
#  (1) HR values changed to 1.5, 1.8, 2.0, 2.5 
#  (2) rho2 values changed to 0.0, 0.1, 0.3, 0.5
#  (3) consider only the combined data (PHS+HPFS) and combined deaths
#      (PCa death+mets)
#  (4) since we will get quintiles based on cases only, the proportion
#      of subjects in the 5-th quintile is 1/5=0.2
#  (5) projected PCa deaths+mets is 310 for all 2412 subjects with blood info
#  (6) projected PCa deaths+mets is 209 for all 1630 subjects with tissue info
#  (7) genotype data can be obtained for subjects with blood info. Hence
#      projected PCa deaths+mets is 310 for all 2412 subjects with genotype info
# 
#
# power for a given sample size
powerEpi.default<-function(n, theta, p, psi, rho2, alpha=0.05)
{
  ua2<-qnorm(1-alpha/2)
  part2<-n*(log(theta))^2*p*(1-p)*psi*(1-rho2)

  power<-pnorm(sqrt(part2)-ua2)

  return(power)
}

# estimate p and rho2 from a pilot data
powerEpi<-function(X1, X2, failureFlag, n, theta, alpha=0.05)
{
  psi<-sum(failureFlag==1)/length(failureFlag)
  p<-sum(X1==1)/length(X1)
  rho<-cor(X1, X2)
  rho2<-rho^2
  
  ua2<-qnorm(1-alpha/2)
  part2<-n*(log(theta))^2*p*(1-p)*psi*(1-rho2)

  power<-pnorm(sqrt(part2)-ua2)

  res<-list(power=power, p=p, rho2=rho2, psi=psi)
  return(res)
}

# sample size required
ssizeEpi.default<-function(power, theta, p, psi, rho2, alpha=0.05)
{
  ua2<-qnorm(1-alpha/2)
  beta<-1-power
  ub<-qnorm(1-beta)

  numer<-(ua2+ub)^2
  denom<-(log(theta))^2*p*(1-p)*psi*(1-rho2)

  n<-ceiling(numer/denom)

  return(n)
}

# estimate p and rho2 from a pilot data
ssizeEpi<-function(X1, X2, failureFlag, power, theta, alpha=0.05)
{
  psi<-sum(failureFlag==1)/length(failureFlag)
  p<-sum(X1==1)/length(X1)
  rho<-cor(X1, X2)
  rho2<-rho^2

  ua2<-qnorm(1-alpha/2)
  beta<-1-power
  ub<-qnorm(1-beta)

  numer<-(ua2+ub)^2
  denom<-(log(theta))^2*p*(1-p)*psi*(1-rho2)

  n<-ceiling(numer/denom)

  res<-list(n=n, p=p, rho2=rho2, psi=psi)
  return(res)
}


# numbe of deaths required
numDEpi.default<-function(power, theta, p, rho2, alpha=0.05)
{
  ua2<-qnorm(1-alpha/2)
  beta<-1-power
  ub<-qnorm(1-beta)

  numer<-(ua2+ub)^2
  denom<-(log(theta))^2*p*(1-p)*(1-rho2)

  D<-ceiling(numer/denom)

  return(D)
}

# numbe of deaths required
numDEpi<-function(X1, X2, power, theta, alpha=0.05)
{

  p<-sum(X1==1)/length(X1)
  rho<-cor(X1, X2)

  ua2<-qnorm(1-alpha/2)
  beta<-1-power
  ub<-qnorm(1-beta)

  numer<-(ua2+ub)^2
  denom<-(log(theta))^2*p*(1-p)*(1-rho^2)

  D<-ceiling(numer/denom)

  res<-list(D=D, p=p, rho2=rho^2)
  return(res)
}



# power calculation for continuous-type variable
powerEpiCont.default<-function(n, theta, sigma2, psi, rho2, alpha=0.05)
{
  ua2<-qnorm(1-alpha/2)
  part2<-n*(log(theta))^2*sigma2*psi*(1-rho2)

  power<-pnorm(sqrt(part2)-ua2)

  return(power)
}

ssizeEpiCont.default<-function(power, theta, sigma2, psi, rho2, alpha=0.05)
{
  ua<-qnorm(1-alpha/2)
  ub<-qnorm(power)

  numer<-(ua+ub)^2
  denom<-(log(theta))^2*sigma2*psi*(1-rho2)

  n<-ceiling(numer/denom)

  return(n)
}


# estimate rho2 based on a pilot data
powerEpiCont<-function(formula, dat, X1, failureFlag, n, theta, alpha=0.05)
{
  sigma2<-var(X1)
  psi<-sum(failureFlag==1)/length(failureFlag)

  res.lm<-lm(formula=formula, data=dat)
  tt<-summary(res.lm)
  rho2<-tt$r.squared

  ua2<-qnorm(1-alpha/2)
  part2<-n*(log(theta))^2*sigma2*psi*(1-rho2)

  power<-pnorm(sqrt(part2)-ua2)

  res<-list(power=power, rho2=rho2, sigma2=sigma2, psi=psi)
  return(res)
}

# estimate rho2 based on a pilot data
ssizeEpiCont<-function(formula, dat, X1, failureFlag, power, theta, alpha=0.05)
{

  sigma2<-var(X1)
  psi<-sum(failureFlag==1)/length(failureFlag)

  res.lm<-lm(formula=formula, data=dat)
  tt<-summary(res.lm)
  rho2<-tt$r.squared


  ua<-qnorm(1-alpha/2)
  ub<-qnorm(power)

  numer<-(ua+ub)^2
  denom<-(log(theta))^2*sigma2*psi*(1-rho2)

  n<-ceiling(numer/denom)

  res<-list(n=n, rho2=rho2, sigma2=sigma2, psi=psi)
  return(res)
}




# p00=Pr(Z1=0, Z2=0)
# p01=Pr(Z1=0, Z2=1)
# p10=Pr(Z1=1, Z2=0)
# p11=Pr(Z1=1, Z2=1)
# n - total number of subjects
# theta = hazard ratio
# psi - proportion of subjects died
# a=n*p00  = no. of Z1=0 and Z2=0
# b=n*p01  = no. of Z1=0 and Z2=1
# c=n*p10  = no. of Z1=1 and Z2=0
# d=n*p11  = no. of Z1=1 and Z2=1
powerEpiInt.default1<-function(n, theta, psi, p00, p01, p10, p11, alpha=0.05)
{
  delta<-1/p00+1/p01+1/p10+1/p11
  ua<-qnorm(1-alpha/2)

  tmp<-sqrt(n*psi*(log(theta))^2/delta)-ua

  power=pnorm(tmp)
  return(power)
}

ssizeEpiInt.default1<-function(power, theta, psi, p00, p01, p10, p11, alpha=0.05)
{
  delta<-1/p00+1/p01+1/p10+1/p11
  ua<-qnorm(1-alpha/2)
  ub<-qnorm(power)

  numer<-(ua+ub)^2*delta
  denom<-log(theta)^2*psi

  n=ceiling(numer/denom)

  return(n)
}

# estimate mya, myb, myc, myd from a pilot study
powerEpiInt<-function(X1, X2, failureFlag, n, theta, alpha=0.05)
{

  psi<-sum(failureFlag==1)/length(failureFlag)
  mya<-sum(X1==0 & X2==0)
  myb<-sum(X1==0 & X2==1)
  myc<-sum(X1==1 & X2==0)
  myd<-sum(X1==1 & X2==1)

  res<-powerEpiInt2(n, theta, psi, 
    mya, myb, myc, myd, alpha)

  res$mya=mya
  res$myb=myb
  res$myc=myc
  res$myd=myd
  res$psi=psi

  return(res)
}



# p00=Pr(Z1=0, Z2=0)
# p01=Pr(Z1=0, Z2=1)
# p10=Pr(Z1=1, Z2=0)
# p11=Pr(Z1=1, Z2=1)
# nobs =mya+myb+myc+myd - total number of subjects
# mya=nobs*p00 - from pilot study
# myb=nobs*p01 - from pilot study
# myc=nobs*p10 - from pilot study
# myd=nobs*p11 - from pilot study
# theta = hazard ratio
# psi - proportion of subjects died
powerEpiInt2<-function(n, theta, psi, 
  mya, myb, myc, myd, alpha=0.05)
{
  ua<-qnorm(1-alpha/2)

  nobs=mya+myb+myc+myd
  # p=Pr(Z1=1) - prevalence of Z1=1
  p=(myc+myd)/nobs

  # q=Pr(Z2=1) - prevalence of Z2=1
  q=(myb+myd)/nobs

  # p0=prevalence of Z1 (Z1=1) in struatum Z2=0
  # p0=Pr(Z1=1 | Z2=0)
  p0<-myc/(mya+myc)

  # p1=prevalence of Z1 (Z1=1) in struatum Z2=1
  # p1=Pr(Z1=1 | Z2=1)
  p1<-myd/(myb+myd)

  numer<-(1-q)*(1-p0)*p0+q*(1-p1)*p1
  numer<-numer^2
  denom<-(1-q)*q*(1-p0)*p0*(1-p1)*p1
  G<-numer/denom

  rho<-(p1-p0)*sqrt(q*(1-q)/(p*(1-p)))
  rho2<-rho^2
  tmp<-sqrt(n*psi*(log(theta))^2*(1-p)*p*(1-rho2)/G)-ua

  power=pnorm(tmp)

  res<-list(power=power, p=p, q=q, p0=p0, p1=p1, rho2=rho2, G=G) 
  return(res)
}

# p and G are estimated based on the data
# mya, myb, myc, myd
ssizeEpiInt2<-function(power, theta, psi, 
  mya, myb, myc, myd, alpha=0.05)
{
  ua<-qnorm(1-alpha/2)
  ub<-qnorm(power)


  nobs=mya+myb+myc+myd
  # p=Pr(Z1=1) - prevalence of Z1=1
  p=(myc+myd)/nobs

  # q=Pr(Z2=1) - prevalence of Z2=1
  q=(myb+myd)/nobs

  # p0=prevalence of Z1 (Z1=1) in struatum Z2=0
  # p0=Pr(Z1=1 | Z2=0)
  p0<-myc/(mya+myc)

  # p1=prevalence of Z1 (Z1=1) in struatum Z2=1
  # p1=Pr(Z1=1 | Z2=1)
  p1<-myd/(myb+myd)

  numer<-(1-q)*(1-p0)*p0+q*(1-p1)*p1
  numer<-numer^2
  denom<-(1-q)*q*(1-p0)*p0*(1-p1)*p1
  G<-numer/denom

  rho<-(p1-p0)*sqrt(q*(1-q)/(p*(1-p)))
  rho2<-rho^2

  numer<-(ua+ub)^2
  denom<-log(theta)^2*psi*(1-p)*p*(1-rho2)

  n<-ceiling(numer*G/denom)

  res<-list(n=n, p=p, q=q, p0=p0, p1=p1, rho2=rho2, G=G) 
  return(res)
}

ssizeEpiInt<-function(X1, X2, failureFlag, power, theta, alpha=0.05)
{

  psi<-sum(failureFlag==1)/length(failureFlag)
  mya<-sum(X1==0 & X2==0)
  myb<-sum(X1==0 & X2==1)
  myc<-sum(X1==1 & X2==0)
  myd<-sum(X1==1 & X2==1)

  res<-ssizeEpiInt2(power, theta, psi, 
  mya, myb, myc, myd, alpha)

  res$mya=mya
  res$myb=myb
  res$myc=myc
  res$myd=myd
  res$psi=psi

  return(res)

}
 
# p=Pr(Z1=1)
powerEpiInt.default0<-function(n, theta, p, psi, G, rho2, 
  alpha=0.05)
{
  ua<-qnorm(1-alpha/2)

  tmp<-sqrt(n*psi*(log(theta))^2*(1-p)*p*(1-rho2)/G)-ua

  power=pnorm(tmp)
  return(power)
}

ssizeEpiInt.default0<-function(power, theta, p, psi, G, rho2, 
  alpha=0.05)
{
  ua<-qnorm(1-alpha/2)
  ub<-qnorm(power)

  numer<-(ua+ub)^2
  denom<-log(theta)^2*psi*(1-p)*p*(1-rho2)

  n<-ceiling(numer*G/denom)

  return(n)
}


