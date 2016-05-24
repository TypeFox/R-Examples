#Convertible Bond price distribution
#t has to be an integer
############################################

setGeneric(
  name="ConvBond_PriceDistribution",
  def = function(.Object,t,T,nCoupons,couponsRate)
  {
    standardGeneric("ConvBond_PriceDistribution")
  }
)
setMethod(
  f="ConvBond_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object,t,T,nCoupons,couponsRate)
  {
    if (!is.numeric(t)||!is.numeric(T)||!is.numeric(nCoupons)||!is.numeric(couponsRate)){stop("All parameters have to be numeric")}
    
    if (floor(t) != t){stop("Parameter t (starting month) must be an integer")}
    
    if (t == T){return (1)}
    
    if (t<=0) 
    { 
      nScenarios <- 1
    }
    else
    {
      nScenarios <- .Object@ParamsScenarios@nScenarios
    }
    
    k <- .Object@ParamsScenarios@k
    vol <-.Object@ParamsScenarios@vol
    
    y <- rep(0,nScenarios)
    
    if (t>=0 && t<=T)
    {
      y <- Bond_PriceDistribution(.Object,t,T,nCoupons,couponsRate)+EuroCall_Stock_PriceDistribution(.Object,t,T,1)
      return(y)
    } else stop("t (starting month) must be positive and lower than maturity T")        
    
    
  }
)