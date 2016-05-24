#Call/Put price distribution on Stock
############################################

setGeneric(
  name="EuroCall_Stock_PriceDistribution",
  def = function(.Object,t,T,Strike)
  {
    standardGeneric("EuroCall_Stock_PriceDistribution")
  }
)
setMethod(
  f="EuroCall_Stock_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object,t,T,Strike)
  { 
    if (!is.numeric(t)||!is.numeric(T)||!is.numeric(Strike)){stop("Parameters t (starting month), T (maturity) and Strike have to be numeric")}
    
    if (floor(t) != t){stop("Parameter t (starting month) must be an integer")}
    
    if (t == T){return (1)}
    
    nn <- length(.Object@stockPaths[1,])
    if (t >= nn || T>=nn) {stop("Not enough yield curve maturities found for this calculation. For the given yield curve, both starting month t=",t," and maturity T=",T," must be <=",round(nn-1,5))}
    
    if (t>0)
    {
        n <- .Object@ParamsScenarios@nScenarios
        St <- .Object@stockPaths[,t]
    }
    else 
    {
        n <-1
        St <- .Object@ParamsScenarios@stock0       
    }
    
    vol <- .Object@ParamsScenarios@vol
    volStock <- .Object@ParamsScenarios@volStock
    rho <- .Object@ParamsScenarios@rho
    k <- .Object@ParamsScenarios@k    
    
    y <- rep(0,n)
    
    if (t>=0 && t<T)
    {
      PtT <- ZCBond_PriceDistribution(.Object,t,T)
      tauT <- Tau(T,vol,volStock,k,rho)
      d1 <- (log(St/(Strike*PtT))+tauT/2)/sqrt(tauT)
      d2 <- d1-sqrt(tauT)
      
      y <- St*pnorm(d1)-Strike*PtT*pnorm(d2)
    } else {stop("t (starting month) must be positive and lower than maturity T")}
    return(y)
  }
)
############################################


setGeneric(
  name="EuroPut_Stock_PriceDistribution",
  def = function(.Object,t,T,Strike)
  {
    standardGeneric("EuroPut_Stock_PriceDistribution")
  }
)
setMethod(
  f="EuroPut_Stock_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object,t,T,Strike)
  { 
    if (!is.numeric(t)||!is.numeric(T)||!is.numeric(Strike)){stop("Parameters t (starting month), T (maturity) and Strike have to be numeric")}
    
    if (floor(t) != t){stop("Parameter t (starting month) must be an integer")}
    
    if (t == T){return (1)}
    
    nn <- length(.Object@stockPaths[1,])
    if (t >= nn || T>=nn) {stop("Not enough yield curve maturities found for this calculation. For the given yield curve, both starting month t=",t," and maturity T=",T," must be <=",round(nn-1,5))}
    
    if (t>0)
    {
      n <- .Object@ParamsScenarios@nScenarios
      St <- .Object@stockPaths[,t]
    }
    else 
    {
      n <-1
      St <- .Object@ParamsScenarios@stock0       
    }
    
    vol <- .Object@ParamsScenarios@vol
    volStock <- .Object@ParamsScenarios@volStock
    rho <- .Object@ParamsScenarios@rho
    k <- .Object@ParamsScenarios@k
    
    y <- rep(0,n)
    
    if (t>=0 && t<T)
    {
      PtT <- ZCBond_PriceDistribution(.Object,t,T)
      tauT <- Tau(T,vol,volStock,k,rho)
      d1 <- (log(St/(Strike*PtT))+tauT/2)/sqrt(tauT)
      d2 <- d1-sqrt(tauT)
      
      y <- -St*pnorm(-d1)+Strike*PtT*pnorm(-d2)
    }  else {stop("t (starting month) must be positive and lower than maturity T")}
    return(y)
  }
)