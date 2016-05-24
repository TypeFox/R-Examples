# Name   : est.R0.SB
# Desc   : Estimation of Reproduction Number distribution using previous R(t)
#          distribution as a Bayesian prior, as presented by Bettencourt & Ribeiro
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

est.R0.SB <- function#Estimate the time dependent reproduction number using a Bayesian approach
### Estimate the time dependent reproduction number using a Bayesian approach. All known data are used as a prior for next iteration (see Details).
##details<< For internal use. Called by est.R0.
##note<< This is the implementation of the method provided by Bettencourt & Ribeiro (2008).
##references<< Bettencourt, L.M.A., and R.M. Ribeiro. "Real Time Bayesian Estimation of the Epidemic Potential of Emerging Infectious Diseases." PLoS One 3, no. 5 (2008): e2185.

(epid, ##<< the epidemic curve
  GT, ##<< generation time distribution
  t=NULL, ##<< Time at which epidemic was observed
  begin=NULL, ##<< At what time estimation begins. Just there for "plot" purposes, not actually used
  end = NULL, ##<< At what time estimation ends. Just there for "plot" purposes, not actually used
  date.first.obs=NULL, ##<< Optional date of first observation, if t not specified
  time.step=1, ##<< Optional. If date of first observation is specified, number of day between each incidence observation
  force.prior = FALSE, ##<< Set to any custom value to force the initial prior as a uniform distribution on [0;value]
  checked=FALSE, ##<< Internal flag used to check whether integrity checks were ran or not.
  ... ##<< parameters passed to inner functions
)
  

# Code

{
  ##details<< Initial prior is an unbiased uniform distribution for R, between 0 and the maximum of incid(t+1) - incid(t). 
  ## For each subsequent iteration, a new distribution is computed for R, using the previous output as new prior.

  #We should think of imports and how to deal with this.
	DNAME =  deparse(substitute(epid))
	CALL = match.call()
  
  #Various class and integrity checks
	if (checked == FALSE) {
	  parameters <- integrity.checks(epid, t, GT, begin, end, date.first.obs, time.step, AR=NULL, S0=NULL, methods="SB")
	  begin <- parameters$begin
	  end <- parameters$end
	}
  epid = check.incid(epid, t, date.first.obs, time.step)
  epid.orig = epid
  begin.nb = which(epid$t == begin)
  end.nb = which(epid$t == end)
  epid$incid <- epid$incid[begin.nb:end.nb]
  epid$t <- epid$t[begin.nb:end.nb]
  time.step <- as.numeric(epid$t[2]-epid$t[1])
  
  ## Beginning of estimation 
	
  Tmax = length(epid$incid)
  gamma = 1/GT$mean # 1/gamma is the infectious period as defined in a standard SIR model
  proba.Rt = list() #list to store successive densities.
  
  
  ## Initial distribution of Rt
  
  if (force.prior == FALSE) {
    #For the first iteration, we admit Rt to be following a unfiorm distribution
    #on [0 ; maximum.delta.incid]
    delta.incid = epid$incid[2:length(epid$incid)] - epid$incid[1:(length(epid$incid)-1)]
    maximum.delta.incid = which.max(delta.incid)
    maximum.delta.incid = epid$incid[maximum.delta.incid+1] - epid$incid[maximum.delta.incid]
  }
  else {
    #Prior by Bettencourt & Ribeiro is an unbiased uniform distribution on [0;3]
    #This one is added for comparison purpose.
    maximum.delta.incid = as.numeric(force.prior)
  }
  
  proba.Rt[[1]] = rep(1/(maximum.delta.incid*100 + 1), (maximum.delta.incid*100 + 1))
  Rt.range = seq(0,maximum.delta.incid, .01)
  most.likely.Rt = rep(0, length(epid$incid))
  
  #Useful results will be stored in a data.frame for easy preview
  Rt.quant = matrix(NA, length(epid$incid), ncol=4)
  colnames(Rt.quant)=c("Date","R(t)", "CI[lower]", "CI[upper]")
  
	
  ## Loop on epidemic duration
	for (s in 2:Tmax) {
    sequence.factorielle <- 1
    if (epid$incid[s] != 0) {
      sequence.factorielle <- seq(from=1, to=epid$incid[s], by=1)
    }
    tmp <- epid$incid[s]*(log(epid$incid[s-1]) + time.step*gamma*(Rt.range - 1)) - epid$incid[s-1]*exp(time.step*gamma*(Rt.range - 1)) - sum(log(sequence.factorielle))
    
    tmp.proba.Rt <- exp(log(proba.Rt[[s-1]]) + tmp)
    #Break if incidence data starts producing NaN values
    if ((sum(tmp.proba.Rt) == 0) | NaN %in% (sum(tmp.proba.Rt))) {
      end.nb <- s-1
      end <- epid$t[end.nb]
      break
    }
    tmp.proba.Rt <- tmp.proba.Rt/sum(tmp.proba.Rt)
    
    proba.Rt[[s]] <- tmp.proba.Rt
    
    # Extracting most likely Rt value.
    most.likely.Rt[s] = Rt.range[which.max(proba.Rt[[s]])]
    
    #Computing Confidence Interval around most plausible value
    ##details<< CI is achieved by a cumulated sum of the R posterior distribution, and corresponds to the 2.5% and 97.5% thresholds
    Rt.quant[s-1,1] <- epid$t[s-1]
    Rt.quant[s-1,2] <- most.likely.Rt[s]
    cumsum.proba <- cumsum(proba.Rt[[s]])
    Rt.quant[s-1,3] <- Rt.range[which.max((cumsum(proba.Rt[[s]])) >= 0.025)]
    Rt.quant[s-1,4] <- Rt.range[which.max((cumsum(proba.Rt[[s]])) >= 0.975)]
  
  }
  
  #Last line of Rt.quant is set to 0 for graphical purposes
  #Rt.quant[length(epid$t),1] = epid$t[length(epid$t)]
  #Rt.quant[length(epid$t),2] = Rt.quant[length(epid$t)-1,2]
  #Rt.quant[length(epid$t),3] = Rt.quant[length(epid$t)-1,3]
  #Rt.quant[length(epid$t),4] = Rt.quant[length(epid$t)-1,4]
  Rt.quant<-data.frame(na.omit(Rt.quant))
  
  #Fix for incorrect storage format in matrix
  if (!is.numeric(epid$t)) {
    Rt.quant$Date <- as.Date(Rt.quant$Date, origin=(as.Date("1970-01-01")))
  }
  
  #Predictive curve is also reported from successive posterior distributions
  pred <- c(epid$incid[1], rep(NA,(end.nb-begin.nb)))
  
  #Based on the most-likely R value at the end of the simulation, will compute the predicted exponential growth curve
  for (t in 1:(Tmax-1)) {
    pred[t+1] <- pred[t]*exp(gamma*time.step*(Rt.quant$R.t[t]-1))
  }
  
  R = Rt.quant[1:(end.nb-begin.nb),2]
	conf.int=Rt.quant[1:(end.nb-begin.nb),3:4]
  rownames(conf.int) = as.character(epid$t[1:(end.nb-begin.nb)])
  return(structure(list(R=R, conf.int=conf.int, proba.Rt=proba.Rt, GT=GT, epid=epid.orig, begin=begin, begin.nb=begin.nb, end=end, end.nb=end.nb, pred=pred, data.name=DNAME, call=CALL, method="Sequential Bayesian", method.code="SB"),class="R0.R"))    #return everything

  ### A list with components:
	### \item{R}{vector of R values.}
	### \item{conf.int}{95% confidence interval for estimates.}
  ### \item{proba.Rt}{A list with successive distribution for R throughout the outbreak.}
	### \item{GT}{Generation time distribution used in the computation.}
	### \item{epid}{Original epidemic data.}
	### \item{begin}{Begin date for the fit.}
  ### \item{begin.nb}{Index of begin date for the fit.}
	### \item{end}{End date for the fit.}
  ### \item{end.nb}{Index of end date for the fit.}
	### \item{pred}{Predictive curve based on most-likely R value.}
  ### \item{data.name}{Name of the data used in the fit.}
  ### \item{call}{Complete call used to generate results.}
	### \item{method}{Method for estimation.}
	### \item{method.code}{Internal code used to designate method.}
}
