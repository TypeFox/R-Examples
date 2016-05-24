# Name   : smooth.Rt
# Desc   : Allows for regrouping real-time reproduction number by a given time period
# Date   : 2011/16/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

smooth.Rt <- function#Smooth real-time reproduction number over larger time period
### Smooth real-time reproduction number over larger time period

(res, ##<< An object of class "R0.R", created by any real-time method (currently implemented: TD and SB)
time.period ##<< Time period to be used for computations.
)
  
  
# Code
  
{
  ##details<<Regrouping Time-Dependant R(t) values, or even Real Time Bayesian most-likely R values (according to R distributions)
  ## should take into account the Generation Time.
  ## Results can be plotted exactly the same was as input estimations, except they won't show any goodness of fit curve.
  
  if (class(res) != "R0.R") {
    stop("Currently, sensitivity analysis from a result object only supports 'R0.R' class objects. Try using res$estimates$TD or res$estimates$SB if they are defined.")
  }
  else if ((res$method %in% c("Time-Dependent","Sequential Bayesian")) == FALSE) {
    stop("Sensitivity analysis can only be conducted on objects with method EG or ML.")
  }
  
  if ((!is.numeric(time.period)) & (!is.integer(time.period))) {
    stop("Error: time.period should be of numeric or integer class.")
  }
  
  #How many groups ?
  nb.groups <- floor(length(res$epid$incid[res$begin.nb:res$end.nb])/time.period)
  
  #New data will have the exact same layout as input
  epid = list(incid = rep(0, nb.groups), t = rep(NA, nb.groups))
  
  Rt.quant <- matrix(NA, nb.groups, ncol=4)
  colnames(Rt.quant) <- c("Date","R(t)", "CI[lower]", "CI[upper]")
  
  #"Date" value is arbitraty set at the first day of the time period considered
  for (t in 1:nb.groups) {
    Rt.quant[t,1] <- res$epid$t[((t-1)*time.period + 1)]
    epid$incid[t] <- sum(res$epid$incid[((t-1)*time.period + 1):(time.period*t)], na.rm=TRUE)
    epid$t[t] <- res$epid$t[((t-1)*time.period + 1)]
    Rt.quant[t,2] <- sum(res$R[((t-1)*time.period + 1):(time.period*t)]*res$epid$incid[((t-1)*time.period + 1):(time.period*t)], na.rm=TRUE)/(sum(res$epid$incid[((t-1)*time.period + 1):(time.period*t)], na.rm=TRUE))
    Rt.quant[t,3] <- sum(res$conf.int[((t-1)*time.period + 1):(time.period*t),1]*res$epid$incid[((t-1)*time.period + 1):(time.period*t)], na.rm=TRUE)/(sum(res$epid$incid[((t-1)*time.period + 1):(time.period*t)], na.rm=TRUE))
    Rt.quant[t,4] <- sum(res$conf.int[((t-1)*time.period + 1):(time.period*t),2]*res$epid$incid[((t-1)*time.period + 1):(time.period*t)], na.rm=TRUE)/(sum(res$epid$incid[((t-1)*time.period + 1):(time.period*t)], na.rm=TRUE))
  }
  
  Rt.quant <- data.frame(Rt.quant)
  
  #Quick fix for imcorrect date format
  if (!is.numeric(res$epid$t)) {
    Rt.quant[,1] <- as.Date(Rt.quant[,1], origin="1970-01-01")
    epid$t <- as.Date(epid$t, origin="1970-01-01")
  }
  
  #Extract required informations for standard return
  R <- Rt.quant[,2]
  conf.int = matrix(data=NA, nrow=dim(Rt.quant)[1], ncol=2)
  colnames(conf.int)=c("lower", "upper")
  conf.int[,1] <- Rt.quant[,3]
  conf.int[,2] <- Rt.quant[,4]
  conf.int<-data.frame(na.omit(conf.int))
  rownames(conf.int) = as.character(Rt.quant[,1])

  
  #return everything
  return(structure(list(R=R, conf.int=conf.int, GT=res$GT, epid=epid, begin=res$begin, begin.nb=1, end=res$end, end.nb=length(epid$t), data.name=res$data.name, call=res$call, method=res$method, method.code=res$method.code),class="R0.R"))

  ### A list with components:
  ### \item{R}{The estimate of the reproduction ratio.}
  ### \item{conf.int}{The 95% confidence interval for the R estimate.}
  ### \item{GT}{Generation time distribution uised in the computation.}
  ### \item{epid}{Original or augmented epidemic data, depending whether impute.values is set to FALSE or TRUE.}
  ### \item{begin}{Starting date for the fit.}
  ### \item{begin.nb}{The number of the first day used in the fit.}
  ### \item{end}{The end date for the fit.}
  ### \item{end.nb}{The number of the las day used for the fit.}
  ### \item{data.name}{The name of the dataset used.}
  ### \item{call}{Call used for the function.}
  ### \item{method}{Method used for fitting.}
  ### \item{method.code}{Internal code used to designate method.}
}
