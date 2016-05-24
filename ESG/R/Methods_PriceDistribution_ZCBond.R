#ZC Bond distribution
#t has to be an integer (starting month)
#T (year)
############################################


setGeneric(
  name="ZCBond_PriceDistribution",
  def = function(.Object,t,T)
  {
    standardGeneric("ZCBond_PriceDistribution")
  }
)
setMethod(
  f="ZCBond_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object,t,T)
  { 
    # t must be an integer (starting month)
    
    if (!is.numeric(t)||!is.numeric(T)){stop("Parameters t (starting month) and T (maturity) have to be numeric")}
    
    if (floor(t) != t){stop("Parameter t (starting month) must be an integer")}
    
    if (t == T){return (1)}
    
    nn <- length(.Object@ForwardRates)
    if (nn <= 0){stop("Zero-coupon curve input must have at least one maturity available")}    
    if (t >= nn || T>=nn) {stop("Not enough yield curve maturities found for this calculation. For the given yield curve, both starting month t=",t," and maturity T=",T," must be <=",round(nn-1,5))}  
        
    # Zero-coupon curve extraction 
    ZC <- .Object@ZCRates #ZCRates[1]=0,ZCRates[2]=R(0,1/12),...
    n <- length(ZC)
        
    # If T is not a multiple of 1/12, then the price is linearly interpolated
    int_un_douze <- floor(T)+seq(from = 0, to = 1, by = 1/12)
    if (sum(T != int_un_douze) == length(int_un_douze))
    {
      t1 <- max(int_un_douze[T >= int_un_douze])
      t2 <- min(int_un_douze[T <= int_un_douze]) 
      maxt2 <- 12*t2
      maxT <- n/12
      
      if (maxt2 <= n)
      {
        if (t1 == 0) {P0t1 <- 1} else {P0t1 <- exp(-ZC[12*t1]*t1)}      
        if (t2 == 0) {P0t2 <- 1} else {P0t2 <- exp(-ZC[12*t2]*t2)} 
        P0T <- P0t1 + ((P0t2-P0t1)/(t2-t1))*(T-t1)
      } else 
        {
          msg <- paste("Not enough yield curve maturities found for this calculation. For the given yield curve, maturity T =",T," must be <=",round(maxT,5))
          stop(msg)
        }            
    }else 
      {
        P0T <- exp(-ZC[12*T]*T)
      }
    
    if (t == 0) 
    {
        return(P0T)             
    } 
    else 
    { 
      if (t>0 && t<=T)
      {
        nScenarios <- .Object@ParamsScenarios@nScenarios
        rt <- .Object@shortRatePaths[,t+1] #shortRatePaths[,1]=r_0,shortRatePaths[,2]=r_1,...
        forward <- .Object@ForwardRates
        
        k <- .Object@ParamsScenarios@k
        vol <-.Object@ParamsScenarios@vol
        
        y <- rep(0,nScenarios)        
        
        P0t <- exp(-ZC[t+1]*t)
        y <- P0T/P0t*exp(-K(T-t,k)*L(t,vol,k)/2+K(T-t,k)*(forward[t+1]-rt))
        
        return(y)
      } else stop("Parameter t (starting month) must be positive and lower than maturity T")        
    }        
  }
)