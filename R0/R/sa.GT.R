# Name   : sa.GT
# Desc   : Finds the possible R0 values with supported estimation methods for
#          varying GT parameters
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

sa.GT=function#Sensitivity analysis of reproduction ratio with varying GT distribution
### Sensitivity analysis of reproduction ratio with varying GT distribution.
##details<< By using different Generation Time (GT) distribution, different estimates of reproduction ratio can be analyzed.

(incid, ##<< incident cases
GT.type, ##<< Type of distribution for GT (see GT.R for details)
GT.mean.range, ##<< mean used for all GT distributions throughout the simulation
GT.sd.range, ##<< Range of standard deviation used for GT distributions. Must be provided as a vector.
begin=NULL, ##<< begin date of the estimation of epidemic
end=NULL, ##<< end date of estimation of the epidemic
est.method, ##<< Estimation method used for sensitivity analysis. Requires a method computing a proper R0 value (and not an instantaneous R(t))
t=NULL, ##<< Dates vector to be passed to estimation function
date.first.obs=NULL, ##<< Optional date of first observation, if t not specified
time.step=1, ##<< Optional. If date of first observation is specified, number of day between each incidence observation
... ##<< parameters passed to inner functions
)
  

# Code
  
{
  #Creating list of GTs and results
  list.GT = list()
  res = list()
  
  #Is estimation method supported?
  if ((est.method %in% c("EG", "ML")) == FALSE) {
    stop("Argument 'est.method' should be any of 'EG' or 'ML' to produce results.")
  }
  
  count = 1
  for (i in 1:length(GT.mean.range)) {
    for (j in 1:length(GT.sd.range)) {
      list.GT[[count]]<-generation.time(GT.type, c(GT.mean.range[[i]], GT.sd.range[j]))
      count = count+1
    }
  }
       
  #tmp.epid is only used to keep trace of (incidence,date) pairs
  tmp.epid = check.incid(incid, t, date.first.obs, time.step)
  
  #All results will be stored in a matrix
  s.a=matrix(NA, nrow=length(list.GT), ncol=6)
  
  #Columns are named so that display is easy to read
  colnames(s.a)=c("GT.Type", "GT.Mean", "GT.SD", "R", "CI.lower", "CI.upper")
  
  for(i in 1:length(list.GT)) {
    #Simulation is ran according to the requested method, with each correct begin/end value
    res=estimate.R(incid, list.GT[[i]], begin=begin, end=end, t=t, date.first.obs=date.first.obs, time.step=time.step, methods=est.method,...)
    

    s.a[i,1]=GT.type
    s.a[i,2]=list.GT[[i]]$mean
    s.a[i,3]=list.GT[[i]]$sd
    if (est.method == "EG") {
      s.a[i,4]=res$estimates$EG$R
      s.a[i,5]=res$estimates$EG$conf.int[1]
      s.a[i,6]=res$estimates$EG$conf.int[2]
    }
    else if (est.method == "ML") {
      s.a[i,4]=res$estimates$ML$R
      s.a[i,5]=res$estimates$ML$conf.int[1]
      s.a[i,6]=res$estimates$ML$conf.int[2]
    }

  }
  
  #Plotting results
  
  return(s.a=s.a)
  
  ### A data frame s.a with following components :
  ### \item{$GT.type}{Distribution law for GT.}
  ### \item{$GT.mean}{Range of means used for tested GTs.}
  ### \item{$GT.sd}{Range of standard deviations used for tested GTs.}
  ### \item{$R}{Computed value for Reproduction Number given GT.type, GT.mean and GT.sd.}
  ### \item{$conf.int[1]}{The lower limit of 95% CI for R.}
  ### \item{$conf.int[2]}{The upper limit of 95% CI for R.}
}
