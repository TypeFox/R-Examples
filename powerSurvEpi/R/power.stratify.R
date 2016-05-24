# v3 created on March 28, 2012
#  (1) we have to assume constant hazard ratio to use Formula (4)
#
# v2 created on March 28, 2012
#  (1) allow for non-constant hazard ratio
#  (2) rename treatment group 2 as treatment group 0 (control group)
#
# created on March 28, 2012
# power calculation for stratified log rank test

#
# Journal of Chronic Diseases
#  Volume 38, Issue 9, 1985, Pages 801-809
#
# Consideration of covariates and stratification in sample size determination for survival time studies
# Mari Palta, Saeid B. Amini

# assuming the treatment effect is constant over time
#  (proportional hazards), and that it is the same in each stratum
#
# exponentail survival with parameter lambda.is in treatment i, stratum s;
# uniform patient entry over the first time unit of the study 
#
#
# T -- total study length
# g.s -- proportion of the total sample size in stratum s
# P.s -- the proportion of stratum s, which is in treatment 1
# lambda.is - hazard for the i-th treatment in stratum s 
# m - number of strata
# HR - hazard ratio (treatment vs control)
# 
# lambda0Vec - hazards for control group in strata 1, ..., m
# HR = lambda1Vec / lambda0Vec
ssize.stratify<-function(power, timeUnit, gVec,
  PVec, HR, lambda0Vec, alpha=0.05, verbose=TRUE
  )
{

  lambda1Vec=lambda0Vec*HR
  
  # gVec - a vector of length m, m=number of strata
  m<-length(gVec)

  # calculate Vs
  # formula (3) on page 802 has 2 typos:
  # there should be a minus sign '-' before lambda.1s and lambda.2s
  # in the first exp function containing '(timeUnit-1)' 
  part1 <- 1-(exp(-lambda1Vec*(timeUnit-1))-exp(-lambda1Vec*timeUnit))/lambda1Vec
  part1 <- PVec*part1
  part2 <- 1-(exp(-lambda0Vec*(timeUnit-1))-exp(-lambda0Vec*timeUnit))/lambda0Vec
  part2 <- (1-PVec)*part2 
  Vs <- part1 + part2

  mu <- log(HR)*sqrt(sum(gVec*PVec*(1-PVec)*Vs))

  za<-qnorm(1-alpha)
  zb<-qnorm(power)
  n = (za+zb)^2/mu^2
  n2 = ceiling(n)
  if(verbose)
  {
    cat("\n*******************\n")
    cat("Number of strata: m=", m, "\n")
    cat("Vs>>\n")
    print(Vs)
    cat("\n")
    cat("mu=", mu, "\n")
    cat("za=", za, ", zb=", zb, "\n")
    cat("n=", n, ", ceiling(n)=", n2, "\n")
    cat("\n*******************\n")
  }
  invisible(n2)
}

tt.power.stratify<-function(power.ini, n, timeUnit, 
gVec,
  PVec, HR, lambda0Vec, alpha=0.05)
{
  n.ini<-ssize.stratify(power=power.ini, timeUnit=timeUnit, 
    gVec=gVec, PVec=PVec, HR=HR, 
    lambda0Vec=lambda0Vec, alpha=alpha, verbose=FALSE
  )
  res<-abs(n.ini-n)
  return(res)
}

power.stratify<-function(n, timeUnit, gVec,
  PVec, HR, lambda0Vec, 
  power.ini=0.8, power.low=0.001, power.upp=0.999, 
  alpha=0.05, verbose=TRUE
  )
{
  res.optim<-stats::optim(par=power.ini, fn=tt.power.stratify, 
      n=n, timeUnit=timeUnit, gVec=gVec,
      PVec=PVec, HR=HR, lambda0Vec=lambda0Vec, alpha=alpha, 
      method = "L-BFGS-B",
      lower = power.low, upper = power.upp)

  if(verbose)
  {
    print(res.optim)
    cat("\npower=", res.optim$par, "\n")
  }

  res<-list(power=res.optim$par, res.optim=res.optim)
  invisible(res)
}

## Assume that a study is to be initiated with 4 years of patient entry
## and 1 year of additional followup. For purposes of the calculations,
##  4 years then constitue one time unit and the total Study length
##  T=1.25 time units (4/4 + 1/4=1.25).
##
## Further assume that patients fall equally into two strata so that
##  g1=g2=0.5 (gVec=c(0.5, 0.5)) and that the 4-year survival probability
##  is 0.10 in stratum 1 and 0.32 in stratum 2. 
## A new treatment is introduced which is expected to increase the survival
## probability to 0.3 in stratum 1.
## With exponential survival, the hazards are 
##  lambda11=-log(0.3)=1.204 (hazard for treatment 1 
##     (new treatment) in stratum 1)
##  and lambda01=-log(0.1)=2.303 (hazard for treatment 2 (traditional treatment)
##     in stratum 1)
##    and the hazard ratio is
##  delta = 2.203/1.204=1.91.
##
## If, in addition, we assume that the treatment effect will be the same in the
## second stratum, the hazards there are 
##  lambda12=1.139/1.91=0.597 (hazard for treatment 1 for stratum 2)
##  and 
##  lambda02=-log(0.32)=1.139 (hazard for treatment 2 for stratum 2)
##  delta = 1.139/0.597=1.91
##
## When randomization is stratified so that P1=P2=0.5, then
##  V1=0.5[1-(exp(-0.25*1.204)-exp(-1.25*1.204))/1.204]
##     +0.5[1-(exp(-0.25*2.303)-exp(-1.25*2.303))]/2.303
##    =0.675
##  V2=0.5[1-(exp(-0.25*1.139)-exp(-1.25*1.139))/1.139]
##     +0.5[1-(exp(-0.25*0.597)-exp(-1.25*0.597))]/0.597
##   =0.451
##
##  mu=log(1.91)*sqrt(0.5*0.675+0.5*0.451) / 2=0.243
##
## alpha=0.05, power=0.9
## n=(1.28+1.64)^2/0.243^2=144
##
#
## type I error rate
#alpha=0.05
#
#power=0.9
#
## assume both strata have hazard ratio 1.91
#delta=1.91
#
## 1.25 time unit (4 year is one time unit)
#T=1.25
#
## assume each of the 2 strata has equal sample size
#gVec=c(0.5, 0.5)
#
## in each stratum, we assume the portion of subjects in treatment 1
##   is the same as that of subjects in treatment 2
#PVec=c(0.5, 0.5)
#
#
#lambda01 = 2.303
#lambda02 = 1.139
#lambda0Vec = c(lambda01, lambda02)
## hazards for treatment 2 in strata 1 and 2
## equal to -log(4-year survival probability) at stratum 1 and stratum 2
## i.e. -log(0.10)=2.303 (stratum 1) and -log(0.32)=1.139 (stratum 2)
## (assuming exponential distribution)
##lambda0Vec=c(2.303, 1.139)
#
## hazards for treatment 1 in strata 1 and 2
##lambda1Vec=c(1.204, 0.597)
#HR<-1/delta # (hazard ratio: treatment vs control)
#cat("HR>>\n"); print(HR); cat("\n");
#lambda1Vec=lambda0Vec/delta
#
#cat("lambda1Vec>>\n"); print(lambda1Vec); cat("\n");
#cat("lambda0Vec>>\n"); print(lambda0Vec); cat("\n");
#
## calculate sample size
#n<-ssize.stratify( 
#  power=power, timeUnit=T, gVec=gVec,
#  PVec=PVec, HR=HR, lambda0Vec=lambda0Vec, 
#  alpha=alpha, verbose=TRUE
#  )
#
#res.power<-power.stratify(n=146, timeUnit=T, gVec=gVec,
#  PVec=PVec, HR=HR, lambda0Vec=lambda0Vec, 
#  power.ini=0.8, power.low=0.001, power.upp=0.999, 
#  alpha=0.05, verbose=TRUE
#  )
#
#
