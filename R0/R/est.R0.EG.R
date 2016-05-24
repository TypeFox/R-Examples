# Name   : est.R0.EG
# Desc   : Estimation of basic Reproduction Number using Exponential Growth method,
#          as presented by Wallinga & Lipsitch
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

est.R0.EG <-function#Estimate R from exponential growth rate
### Estimate R from exponential growth rate.
##note<< This is the implementation of the method provided by Wallinga & Lipsitch (2007).
##details<< For internal use. Called by est.R0.
##references<< Wallinga, J., and M. Lipsitch. "How Generation Intervals Shape the Relationship Between Growth Rates and Reproductive Numbers." Proceedings of the Royal Society B: Biological Sciences 274, no. 1609 (2007): 599.

( epid, ##<< object containing epidemic curve data. see Details.
  GT, ##<< generation time distribution
  t=NULL, ##<< Vector of dates at which incidence was calculated
  begin=NULL, ##<< At what time estimation begins 
  end=NULL, ##<< Time at which to end computation
  date.first.obs=NULL, ##<< Optional date of first observation, if t not specified
  time.step=1, ##<< Optional. If date of first observation is specified, number of day between each incidence observation
  reg.met="poisson", ##<< Regression method used. Default is "poisson" (for GLM), but can be forced to "linear".
  checked=FALSE, ##<< Internal flag used to check whether integrity checks were ran or not.
  ... ##<< parameters passed to inner functions
)


# Code

{
  DNAME =  deparse(substitute(epid))
  #CALL = sapply(match.call()[-1], deparse)
  CALL = match.call()
  
  # Various class and integrity check
  if (checked == FALSE) {
    parameters <- integrity.checks(epid, t, GT, begin, end, date.first.obs, time.step, AR=NULL, S0=NULL, methods="EG")
    begin <- parameters$begin
    end <- parameters$end
  }
  epid = check.incid(epid, t, date.first.obs, time.step)
  begin.nb = which(epid$t == begin)
  end.nb = which(epid$t == end)

  
  
  #Backup original epidemic data
  epid.orig = epid
  
	#Epidemic data used for fit is truncted to keep only values within [begin,end]
  #epid = list(incid=epid$incid[begin.nb:end.nb], t=epid$t[begin.nb:end.nb], t.glm=seq(from=begin.nb, to=end.nb, by=1))
  #epid = list(incid=epid$incid[begin.nb:end.nb], t=epid$t[begin.nb:end.nb], t.glm=epid$t[begin.nb:end.nb])
  epid = list(incid=epid$incid[begin.nb:end.nb], t=epid$t[begin.nb:end.nb])

  # Different methods to estimate epidemic growth rate (r) from data
  # Method 1 == Linear regression

  ##details<< method "poisson" uses Poisson regression of incidence.
  ## method "linear" uses linear regression of log(incidence)
  reg.met = match.arg(reg.met)
  if (reg.met == "linear") {
     #tmp <-lm((log(incid)) ~ t.glm, data=epid)
     tmp <-lm((log(incid)) ~ t, data=epid)
     Rsquared = summary(tmp)$r.squared
     r <- coefficients(tmp)[2]
     conf.int = confint(tmp)[2,]
	pred = exp(predict(tmp,type="response"))
  } 
  
  # Method 2 == Poisson regression
  else if (reg.met == "poisson") {
    #tmp <- glm(incid ~ t.glm, family=poisson(), data=epid)
    tmp <- glm(incid ~ t, family=poisson(), data=epid)
    Rsquared = (tmp$null.deviance-tmp$deviance)/(tmp$null.deviance)
    r <- coefficients(tmp)[2]
    confint = confint(tmp)[2,]
    pred= predict(tmp,type="response")
  }

	# GT represents Generation Time distribution -> #WE DEFINE A GI type SO THAT NO CHECK REQUIRED BUT type/class
	# Apply method of Wallinga / Lipsitch for discretized Laplace transform
	
	R = as.numeric(R.from.r(r,GT))
  ##details<< CI is computed from the 1/M(-r) formula using bounds on r from the Poisson regression.
	R.inf = as.numeric(R.from.r(confint[1],GT))
	R.sup = as.numeric(R.from.r(confint[2],GT))
	
  #return(structure(list(R=R, conf.int = c(R.inf,R.sup), r=r, conf.int.r=confint,epid=epid.orig, GT=GT, data.name =DNAME, begin=begin,end=end,method="Exponential Growth",pred=pred,fit=met),class="R"))
  return(structure(list(R=R, conf.int = c(R.inf,R.sup), r=r, conf.int.r=confint, Rsquared=Rsquared, epid=epid.orig, GT=GT, data.name=DNAME, call=CALL, begin=begin, begin.nb=begin.nb, end=end, end.nb=end.nb, method="Exponential Growth",pred=pred,fit=reg.met, method.code="EG"),class="R0.R"))
  
  ### A list with components:
  ### \item{R}{The estimate of the reproduction ratio.}
  ### \item{conf.int}{The 95% confidence interval for the R estimate.}
  ### \item{r}{Exponential growth rate of the epidemic.}
  ### \item{conf.int.r}{Confidence interval of the exponential growth rate of the epidemic.}
  ### \item{Rsquared}{The deviance R-squared measure for the considered dates and model.}
  ### \item{epid}{Original epidemic data.}
  ### \item{GT}{Generation time distribution uised in the computation.}
  ### \item{data.name}{Name of the data used in the fit.}
  ### \item{begin}{Starting date for the fit.}
  ### \item{begin.nb}{The number of the first day used in the fit.}
  ### \item{end}{The end date for the fit.}
  ### \item{end.nb}{The number of the las day used for the fit.}
  ### \item{fit}{Method used for fitting.}
  ### \item{pred}{Prediction on the period used for the fit.}
  ### \item{method}{Method for estimation.}
  ### \item{method.code}{Internal code used to designate method.}
}
