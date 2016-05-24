# Name   : est.R0.TD
# Desc   : Estimation of time-Dependent Reproduction Number using an infectious-infected
#          network, as presented by Wallinga & Teunis
# Date   : 2011/11/09
# Author : Boelle, Obadia, Cauchemez
###############################################################################


# Function declaration

est.R0.TD <- function#Estimate the time dependent reproduction number
### Estimate the time dependent reproduction number according to Wallinga & Teunis.
##details<< For internal use. Called by est.R0.
##note<< This is the implementation of the method provided by Wallinga & Teunis (2004). Correction for estimation in real time is implemented as in Cauchemez et al, AJE (2006).
##references<< Wallinga, J., and P. Teunis. "Different Epidemic Curves for Severe Acute Respiratory Syndrome Reveal Similar Impacts of Control Measures." American Journal of Epidemiology 160, no. 6 (2004): 509.

(epid, ##<< epidemic curve.
  GT, ##<< generation time distribution.
  import=NULL, ##<< Vector of imported cases.
  n.t0=NULL, ##<< Number of cases at time 0.
  t=NULL, ##<< Vector of dates at which incidence was measured. 
  begin=NULL, ##<< At what time estimation begins. Just here for "plot" purposes, not actually used
  end = NULL, ##<< At what time estimation ends. Just here for "plot" purposes, not actually used
  date.first.obs=NULL, ##<< Optional date of first observation, if t not specified
  time.step=1, ##<< Optional. If date of first observation is specified, number of day between each incidence observation.
  q = c(0.025, 0.975), ##<< Quantiles for R(t). By default, 5% and 95%
  correct=TRUE, ##<< Correction for cases not yet observed (real time).
  nsim = 10000, ##<< Number of simulations to be run to compute quantiles for R(t)
  checked=FALSE, ##<< Internal flag used to check whether integrity checks were ran or not.
  ... ##<< parameters passed to inner functions
)
  

# Code

{

  if (nsim < 1000) warning("Accurate confidence interval for R(t) requires a large number of simulations. Consider increasing 'nsim'")
  if (nsim>100) warning("Simulations may take several minutes.")
  
	DNAME =  deparse(substitute(epid))
	CALL = match.call()
  
  #Various class and integrity checks
  #transforms epidemic to the right format
  if (checked == FALSE) {
    parameters <- integrity.checks(epid, t, GT, begin, end, date.first.obs, time.step, AR=NULL, S0=NULL, methods="TD")
    begin <- parameters$begin
    end <- parameters$end
  }
	epid = check.incid(epid, t, date.first.obs, time.step)
  begin.nb = which(epid$t == begin)
  end.nb = which(epid$t == end)
  epid.bak <- epid
  # en dessous : devrait etre uniquement a la fin 
  #epid$incid <- epid$incid[begin.nb:end.nb]
  #epid$t <- epid$t[begin.nb:end.nb]
  
  # Warning if length(GT) <= longest gap in epidemic curve
  t <- diff(c(FALSE, epid$incid==0, FALSE), 1)
  start <- which(t==1)
  end <- which(t==-1)
  if (length(start) > 0 & length(end) > 0) { 
    longest <- max(end-start)
    if (longest > length(GT$GT)) warning(paste("Gap in epidemic curve is longer than the generation interval. Consider using a different GT distribution (maybe with \"truncate=", longest, "\" (length of longest gap))."), sep="")
  }
  
  #Imported cases should be provided as a vector of the same length as incid.
  #If no imported cases is provided, import is set to 0.
  if (is.null(import)) {
    import<-rep(0, length(epid$incid))
  }
  ##note<<If imported cases are provided, they are counted in addition 
  ## to autonomous cases. The final plot  
  ## will show overall incidence.
  
  if (!is.null(import) & (length(import) != length(epid$incid))) {
    stop("Vector of imported cases should have the same length as 'epid' data.")
  }
  
  #Unknown initial n0 is set to first recorded value
	if (is.null(n.t0)) {
    n.t0=epid$incid[1]
    warning("Using initial incidence as initial number of cases.")
	}
  
  #Initial n0 larger than first recorded value is an error
	if (n.t0 > epid$incid[1])  {
		stop(paste("Provided initial number of cases (n.t0=",n.t0,") is larger than incidence on begin day (=",epid$incid[begin],")"))
  }

  #Beginning of estimation 
	Tmax = length(epid$incid)
	GT.pad = GT$GT
	if (length(GT.pad) < Tmax) {
		#Pad GI with 0s to allow for matrix multiplication
		GT.pad <- c(GT.pad, rep(0,Tmax - length(GT.pad)))
	}
  
	#Matrix to store probabilities
  #P[i,j] the mean number of offspring at time j caused by infected at time i
  #p[i,j] is proba that a case i infected case j
	P <- matrix(0,ncol=Tmax,nrow=Tmax)
	p <- matrix(0,ncol=Tmax,nrow=Tmax)
  
  
  #At the same time, multiple simulation to estimate quantiles for R(t) values.
  #At each time unit 's', we count how many offspring come from each of [1:s] time unit 
  multinom.simu = vector("list", Tmax)
  multinom.simu[[1]] = matrix(0, Tmax, nsim)
  
  if (epid$incid[1]-n.t0 > 0) {
		# non index cases allocated to one of other cases at day t0, including other non index
		P[1,1]<-(epid$incid[1]-n.t0)/(epid$incid[1]-1)
    p[1,1]<-1
		multinom.simu[[1]][1,] = rmultinom(nsim, epid$incid[1]-n.t0, p[1:1,1])
	}
  
  #imported cases are added to autonomous so that they can be "infectors"infectors" after their arrival
  epid.orig<-epid
  epid$incid = epid$incid + import
	
  #Loop on epidemic duration
	for (s in 2:Tmax) {
    
    multinom.simu[[s]] = matrix(0, Tmax, nsim)
    
		#If autonomous cases were incident on day s, we have to find how previous cases may have contributed to them
		if ((epid$incid[s]-import[s]>0)) {
      #weights for index case detected on day s
			weight.cases.for.s <-(epid$incid[1:s]-c(rep(0,s-1),1+import[s]))*GT.pad[s:1]
			
      #normalization
			weight.cases.for.s <- weight.cases.for.s / sum(weight.cases.for.s)
			
      #Likelihood of a infected at time s to have been caused by parent in [1:s]
      #prob.cases.for.s[i] is the expected number of offspring on day s from one individual on day i
      prob.cases.for.s <- weight.cases.for.s * (epid$incid[s]-import[s])/(epid$incid[1:s]-c(rep(0,s-1),1+import[s]))
			prob.cases.for.s[epid$incid[1:(s-1)]==0] <- 0
			
      #Should not be required, but just in case
      #only 1 icident case on day s: can't be its own infector
			if(epid$incid[s]-import[s] == 1) {
        prob.cases.for.s[s] <- 0
			}
      
      #P will store the average number of offspring coming from each time unit :
      #P[i,j] is the mean number of offspring at time j caused by infected at time i
			P[1:s,s]<-prob.cases.for.s
      
      #Probabilities (empirical) are obtained obtained from weight.cases
      p[1:s,s]<-weight.cases.for.s
      
      #We generate a high number of simulations with the computed probabilities
      #to find out how many infection result from each time unit.
      #multinom.sim list is updated from previous value, so we don't have to sum all lists at the end
      multinom.simu[[s]][1:s,] = multinom.simu[[s-1]][1:s,] + rmultinom(nsim,epid$incid[s]-import[s],p[1:s,s])
		}
    else {
      P[1:s,s]<-0
      p[1:s,s]<-0
      multinom.simu[[s]][1:s,] = multinom.simu[[s-1]][1:s,]
    }
	}
  
  
  #We now have enough data to compute R (from infection network P)
  #along with its 5% and 95% quantiles (from multiple simulations and p)
	R.WT<-apply(P,1,sum)    # Wallinga and Teunis definition
	R.corrected<-R.WT/(cumsum(GT.pad[1:Tmax]))[Tmax:1]   # Corrected for real time
  if (is.na(R.corrected[length(epid$incid)])) {
    R.corrected[length(epid$incid)] <- 0
  }
  
  #Simulated incidence at each time unit is the sum of all cases,
  #stored in the last element of multinom.sim list
  total.infected.by.time.unit.simu = multinom.simu[[length(epid$incid)]]
  R.simu<-total.infected.by.time.unit.simu/c(epid$incid)
  R.simu.corrected<-R.simu/(cumsum(GT.pad[1:Tmax]))[Tmax:1] # Corrected for real time
  R.simu.corrected[Tmax,] <- 0
  
  #Initiating quantile matrix
  quant.simu = matrix(0, Tmax, 2)
  quant.simu.corrected = matrix(0, Tmax, 2)
  
  #Quantile are computed with the nsim*length(incid) values of R(t)
  for (s in 1:Tmax) {
    if (epid$incid[s] == 0) {
      R.WT[s] <- 0
      R.simu[s] <- 0
      R.corrected[s] <- 0
      R.simu.corrected[s] <- 0
    }
    quant.simu[s,] = quantile(R.simu[s,], q, na.rm=TRUE)
    quant.simu.corrected[s,] = quantile(R.simu.corrected[s,], q, na.rm=TRUE)
  }
  
  # changed Tmax to end
  ##details<< CI is computed by multinomial simulations at each time step with the expected value of R.
  conf.int = matrix(data=NA, nrow=end.nb, ncol=2)
  colnames(conf.int)=c("lower", "upper")
  
	if (correct == TRUE) {
	  R = R.corrected[begin.nb:end.nb]
    conf.int[begin.nb:end.nb,1] = quant.simu.corrected[begin.nb:end.nb,1]
	  conf.int[begin.nb:end.nb,2] = quant.simu.corrected[begin.nb:end.nb,2]
	}
  
  else {
    R = R.WT[begin.nb:end.nb]
    conf.int[begin.nb:end.nb,1] = quant.simu[begin.nb:end.nb,1]
    conf.int[begin.nb:end.nb,2] = quant.simu[begin.nb:end.nb,2]
	}
  names(R) = epid$t[begin.nb:end.nb]
	
	conf.int<-data.frame(na.omit(conf.int))
  	rownames(conf.int) = as.character(epid$t[begin.nb:end.nb])

	#Fix for incorrect storage format in matrix
#  if (!is.numeric(epid$t)) {
#    Rt.quant$Date <- as.Date(Rt.quant$Date, origin=(as.Date("1970-01-01")))
#  }
  
	#Now we generate the theoretical incidence graph to compare with the model
  pred = epid$incid
  pred[2:(length(epid$incid)+length(GT$GT))] = 0

  for (s in 1:end.nb) {
    pred[s:(s+length(GT$GT)-1)] = pred[s:(s+length(GT$GT)-1)] + R[s] * epid$incid[s] * GT$GT
  }
	pred = pred[1:end.nb]
	
	return(structure(list(R=R,conf.int=conf.int, P=P, p=p, GT=GT, epid=epid.bak, import=import, pred=pred, begin=begin, begin.nb=begin.nb, end=end, end.nb=end.nb, data.name=DNAME, call=CALL, method="Time-Dependent", method.code="TD"),class="R0.R"))    #return everything
	
  ### A list with components:
  ### \item{R}{vector of R values.}
  ### \item{conf.int}{95% confidence interval for estimates.}
	### \item{P}{Matrix of who infected whom.}
  ### \item{p}{Probability of who infected whom (values achieved by normalizing P matrix).}
  ### \item{GT}{Generation time distribution used in the computation.}
	### \item{epid}{Original epidemic data.}
  ### \item{import}{Vector of imported cases.}
  ### \item{pred}{Theoretical epidemic data, computed with estimated values of R.}
  ### \item{begin}{Starting date for the fit.}
  ### \item{begin.nb}{The number of the first day used in the fit.}
  ### \item{end}{The end date for the fit.}
  ### \item{end.nb}{The number of the las day used for the fit.}
  ### \item{data.name}{Name of the data used in the fit.}
  ### \item{call}{Call used for the function.}
	### \item{method}{Method for estimation.}
  ### \item{method.code}{Internal code used to designate method.}
}
