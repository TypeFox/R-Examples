#Bond price distribution
#t has to be an integer
############################################


setGeneric(
  name="Bond_PriceDistribution",
  def = function(.Object,t,T,nCoupons,couponsRate)
  {
    standardGeneric("Bond_PriceDistribution")
  }
)
setMethod(
  f="Bond_PriceDistribution",
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
      y <- ZCBond_PriceDistribution(.Object,t,T) #P(t,T)
      betat <- Beta(nCoupons,t,T)
      
      for (i in 0:(betat-1))
      {
        alphai <- Alpha(nCoupons,i,T) #Coupon date
        y <- y + couponsRate*ZCBond_PriceDistribution(.Object,t,alphai) #c*P(t,alpha_i)
      }
      y <- y + CouponTimeValue(couponsRate,nCoupons,t,T) #Coupon's time value addition
      
      return(y)
    } else stop("t (starting month) must be positive and lower than maturity T")        
    
    
  }
)