# created on Feb. 25, 2009 by Weiliang Qiu
#
# Rosner B. (2006). Fundamentals of Biostatistics (6th edition). Thomson Brooks/Cole.
# Section 14.12 page 807
#

# Estimation of Power for the Comparison of Survival Curves Between Two Groups under
# the Cox Proportional-Hazards Model (for clinical trials)
#
# Suppose we want to compare the survival curves between an experimental group (E) and 
# a control group (C) in a clinical trial with nE participants in the E group
# and nC participants in the C group, with a maximum follow-up of t years. We wish to 
# test the hypothesis H0: RR=1 versus H1: RR not equal to 1, where RR=underlying hazard ratio
# for the E group versus the C group. We postulate a hazard ratio of RR under H1, and want to 
# conduct a two-sided test with significance level alpha. If the ratio of participants in group 
# 1 compared to group 2 = nE/nC=k, then the power of the test is
# power=Phi(sqrt(k*m)*|RR-1|/(k*RR+1)-za), where za is the 100 (1-alpha/2) percentile of
# the standard normal distribution N(0, 1), Phi is the CDF of N(0, 1).
#
# The method was proposed in Freedman, L.S. (1982). Tables of the number of patients required in 
# clinical trials using the log-rank test. Statistics in Medicine, 1, 121-129.

# k=ratio of participants in group 1 compared to group 2
# m=expected total number of events over both groups
# RR=postulated hazard ratio
# alpha=type I error rate
# pC=probability of failure in group C over the maximum time period of the study (t years)
# pE=probability of failure in group E over the maximum time period of the study (t years)

# power calculation for clinical trial studies
powerCT.default<-function(nE, nC, pE, pC, RR, alpha=0.05)
{
  k<-nE/nC
  m<-nE*pE+nC*pC

  res<-powerCT.default0(k, m, RR, alpha)

  return(res)
}


powerCT.default0<-function(k, m, RR, alpha=0.05)
{
  za<-qnorm(1-alpha/2)

  tt<-sqrt(k*m)*abs(RR-1)/(k*RR+1)-za

  power<-pnorm(tt)
  return(power)
}

# sample size calculation
ssizeCT.default<-function(power, k, pE, pC, RR, alpha=0.05)
{
  za<-qnorm(1-alpha/2)
  zb<-qnorm(power)
 
  m=((k*RR+1)/(RR-1))^2*(za+zb)^2/k

  denom<-k*pE+pC
  nE=ceiling(m*k/denom)
  nC=ceiling(m/denom)

  res<-c(nE, nC)
  names(res)<-c("nE", "nC")
  return(res)
}

# get lambda_i and delta_i
# lambdai=approximate hazard at time i in the C group, i=1, ..., t
#        = Pr(failure at time i among participants in the C group | a participant
#             has survived to time i-1 and is ot censored at time i-1)
# deltai=Pr(a participant is cxensored at time i | he was followed up to time i and has not failed),
#       i=1,...,t           
# deltai is assumed the same in each group

# dat - a data frame with at least three columns: time, status, and x
# formula - e.g. Surv(time, status) ~ x, where x indicates which group a subject is in
# x can take only two possible values: "C" (control group) and "E" (Treatment group)
# and x should be a factor object of R
getLambdaDetla<-function(formula, dat, RR=1.5)
{
  fit<-survfit(formula=formula, data=dat)

  # Pr (survive to time t_i | survived up to time t_{i-1})
  surv<-fit$surv
  strata<-fit$strata

  str<-names(strata)
  str1<-unlist(strsplit(str[1], split="="))[2]
  str2<-unlist(strsplit(str[2], split="="))[2]
  if(str1 != "C" || str2 != "E")
  {
    stop("group variable should take only two possible values: C and E")
  }

  if(str1 == "E")
  { nE<-strata[1]
    nC<-strata[2]
    nTotal<-nE+nC
    set<-1:nTotal
    pos.E<-1:nE
    pos.C<-set[-pos.E]
    #surv.E<-surv[1:nE]
    #surv.C<-surv[-(1:nE)]
  } else if (str1 == "C") {
    nE<-strata[2]
    nC<-strata[1]
    nTotal<-nE+nC
    set<-1:nTotal
    pos.C<-1:nC
    pos.E<-set[-pos.C]

    #surv.E<-surv[1:nE]
    #surv.C<-surv[-(1:nE)]
  }

  surv.E<-surv[pos.E]
  surv.C<-surv[pos.C]

  # number of failed at time t
  #n.event.C<-fit$n.event[-c(1:nE)]
  n.event.C<-fit$n.event[pos.C]
  # number of at risk at time t 
  #  = number of failed  at time t  
  #    + number of alived at time t 
  #    + number of censored at time t
  #n.risk.C<-fit$n.risk[-c(1:nE)]
  n.risk.C<-fit$n.risk[pos.C]

  nTimes<-length(surv.C)
  nTimes1<-nTimes-1
  n.censored.C<-rep(0, nTimes)

  for(i in  1:nTimes1)
  { 
    # number of censored + number of  alived
    tmp <-n.risk.C[i]-n.event.C[i] 
    n.censored.C[i] <-tmp-n.risk.C[i+1]
  }
  n.censored.C[nTimes]<-n.risk.C[nTimes]-n.event.C[nTimes]
  n.survive.C<-n.risk.C-n.event.C-n.censored.C

  lambda.C<-rep(0, nTimes)
  lambda.C[1]<- n.survive.C[1]/n.risk.C[1]
  for(i in 2:nTimes)
  {
    lambda.C[i]<-(n.survive.C[i-1]-n.event.C[i])/n.survive.C[i-1]
  }
  lambda.C<-1-lambda.C

  RRlambda.C<-RR*lambda.C

  #times<-fit$time[1:nC]
  times<-fit$time[pos.C]
  mat.event<-cbind(times, n.event.C, n.censored.C, n.survive.C, n.risk.C)
  mat.event<-rbind(c(0, 0, 0, 0, n.risk.C[1]), mat.event)
  colnames(mat.event)<-c("time", "nEvent.C", "nCensored.C", "nSurvive.C", "nRisk.C")


  delta.C<-n.censored.C/(n.censored.C+n.survive.C)


  lambda.C1<-lambda.C[1]
  A<-rep(0, nTimes)
  B<-rep(0, nTimes)
  C<-rep(0, nTimes)

  A[1]<-1
  B[1]<-1
  C[1]<-1
  for(i in 2:nTimes) 
  {
    i1<-i-1 
    A[i]<-prod(1-lambda.C[1:i1])
    B[i]<-prod(1-RRlambda.C[1:i1])
    C[i]<-prod(1-delta.C[1:i1])
  }

  D<-lambda.C*A*C
  E<-RRlambda.C*B*C

  pC<-sum(D)
  pE<-sum(E)
  

  mat<-cbind(times, lambda.C, RRlambda.C, delta.C, A, B, C, D, E)
  mat2<-rbind(c(0, 0, 0, 0, NA, NA, 1.0, NA, NA), mat)
  colnames(mat2)<-c("time", "lambda", "RRlambda", "delta", "A", "B", "C", "D", "E")

  res<-list(mat.lambda=mat2, mat.event=mat.event, pC=pC, pE=pE)

  return(res)
}

# power calculation. Estimates of pE and pC are based on 'dat'
powerCT<-function(formula, dat, nE, nC, RR, alpha=0.05)
{

  res<-getLambdaDetla(formula, dat, RR)
  pE<-res$pE
  pC<-res$pC

  k<-nE/nC
  m<-nE*pE+nC*pC

  power<-powerCT.default0(k, m, RR, alpha)

  res$power<-power
  return(res)
}

# sample size calculation. Estimates of pE and pC are based on 'dat'
ssizeCT<-function(formula, dat, power, k, RR, alpha=0.05)
{

  res<-getLambdaDetla(formula, dat, RR)
  pE<-res$pE
  pC<-res$pC

  ssize<-ssizeCT.default(power, k, pE, pC, RR, alpha)

  res$ssize<-ssize

  return(res)
}


