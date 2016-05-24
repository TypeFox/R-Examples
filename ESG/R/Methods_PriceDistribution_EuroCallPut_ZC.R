#ZC bond european call/put price distribution
#t might be an integer
#s,T have to be multiple of 1/12 (monthly time)
############################################


setGeneric(
  name="EuroCall_ZC_PriceDistribution",
  def = function(.Object,t,T,s,Strike)
  {
    standardGeneric("EuroCall_ZC_PriceDistribution")
  }
)
setMethod(
  f="EuroCall_ZC_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object,t,T,s,Strike)
  { 
    if (!is.numeric(t)||!is.numeric(T)||!is.numeric(s)||!is.numeric(Strike)){stop("Parameters t (starting month), T, s (maturities) and Strike have to be numeric")}
    
    if (floor(t) != t){stop("Parameter t (starting month) must be an integer")}
    
    if (t == 0 && T == 0){return (1)}
    
    nn <- length(.Object@ZCRates)/12
    if (nn <= 0){stop("Zero-coupon curve input must have at least one maturity available")}    
    if (T>=nn || s>=nn) {stop("Not enough yield curve maturities found for this calculation. For the given yield curve, starting month t=",t," and maturites T=",T," and s=",s," must be <=",round(nn-1,5))}    
        
    if (t==0) 
    {
      n <- 1
      y <- 0
    }
    else{
      n <- .Object@ParamsScenarios@nScenarios
      y <- rep(0,n)
    }
    
    vol <- .Object@ParamsScenarios@vol
    k <- .Object@ParamsScenarios@k
    
    
    if (t<T & T<=s)
    {
      H <- sqrt(vol^2/(2*k^3)*((exp(-k*(s-T))-1)^2-(exp(-k*s)-exp(-k*T))^2))
      d1 <- rep(0,n)
      Pts <- ZCBond_PriceDistribution(.Object,t,s)
      PtT <- ZCBond_PriceDistribution(.Object,t,T)
      d1 <- 1/H*log(Pts/(PtT*Strike)+H/2)
      d2 <- d1-H
      
      y <- Pts*pnorm(d1)-Strike*PtT*pnorm(d2)
    } else stop("There must be t < T <= s")
    return(y)
  }
)
############################################


setGeneric(
  name="EuroPut_ZC_PriceDistribution",
  def = function(.Object,t,T,s,Strike)
  {
    standardGeneric("EuroPut_ZC_PriceDistribution")
  }
)
setMethod(
  f="EuroPut_ZC_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object,t,T,s,Strike)
  {
    if (!is.numeric(t)||!is.numeric(T)||!is.numeric(s)||!is.numeric(Strike)){stop("Parameters t (starting month), T, s (maturities) and Strike have to be numeric")}
    
    if (floor(t) != t){stop("Parameter t (starting month) must be an integer")}
    
    if (t == 0 && T == 0){return (1)}
    
    nn <- length(.Object@ZCRates)/12
    if (nn <= 0){stop("Zero-coupon curve input must have at least one maturity available")}    
    if (T>=nn || s>=nn) {stop("Not enough yield curve maturities found for this calculation. For the given yield curve, starting month t=",t," and maturites T=",T," and s=",s," must be <=",round(nn-1,5))}    
        
    if (t==0) 
    {
      n <- 1
      y <- 0
    }
    else{
      n <- .Object@ParamsScenarios@nScenarios
      y <- rep(0,n)
    }
    
   
    
    if (t<T & T<=s)
    {
      vol <- .Object@ParamsScenarios@vol
      k <- .Object@ParamsScenarios@k
      H <- sqrt(vol^2/(2*k^3)*((exp(-k*(s-T))-1)^2-(exp(-k*s)-exp(-k*T))^2))
      d1 <- rep(0,n)
      Pts <- ZCBond_PriceDistribution(.Object,t,s)
      PtT <- ZCBond_PriceDistribution(.Object,t,T)
      d1 <- 1/H*log(Pts/(PtT*Strike)+H/2)
      d2 <- d1-H
      
      y <- -Pts*pnorm(-d1)+Strike*PtT*pnorm(-d2)
    } else stop("There must be t < T <= s")
    return(y)
  }
)