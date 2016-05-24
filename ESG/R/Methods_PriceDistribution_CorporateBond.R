#Corporate bond price distribution
#t integer
#Omega represents the part recovered in case of default
############################################

setGeneric(
  name="CBond_PriceDistribution",
  def = function(.Object,t,T,nCoupons,couponsRate,omega)
  {
    standardGeneric("CBond_PriceDistribution")
  }
)
setMethod(
  f="CBond_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object,t,T,nCoupons,couponsRate,omega)
  {     
    if (!is.numeric(t)||!is.numeric(T)||!is.numeric(nCoupons)||!is.numeric(couponsRate)||!is.numeric(omega)){stop("All the parameters have to be numeric")}
    
    if (floor(t) != t){stop("Parameter t (starting month) must be an integer")}
    
    if (t == T){return (1)}
    
    n <- length(.Object@liquiditySpreadPaths[1,])
    if (t>=n) 
    {
      msg <- paste("Not enough yield curve maturities found for this calculation. For the given yield curve, starting month t=",t," must be <=",round(n-1,5))
      stop(msg)
    }
    
    if (t<=0) 
    { 
      nScenarios <- 1
      y <- 0
      y1 <- 0
      y2 <- 0
      y3 <- 0
      gammat <- .Object@liquiditySpreadPaths[1,t+1]
      lambdat <- .Object@defaultSpreadPaths[1,t+1]
    }
    else
    {
      nScenarios <- .Object@ParamsScenarios@nScenarios
      y <- rep(0,nScenarios)
      y1 <- rep(0,nScenarios)
      y2 <- rep(0,nScenarios)
      y3 <- rep(0,nScenarios)
      gammat <- .Object@liquiditySpreadPaths[,t+1]
      lambdat <- .Object@defaultSpreadPaths[,t+1]
    }    
    
    if (t>=0 && t<=T)
    {
      alpha <- .Object@ParamsScenarios@alpha
      beta <- .Object@ParamsScenarios@beta
      eta <- .Object@ParamsScenarios@eta
      volDefault <- .Object@ParamsScenarios@volDefault
      betat <- Beta(nCoupons,t,T)
        
        #Calcul du premier terme
        for (i in 0:(betat-1))
        {
          alphai <- Alpha(nCoupons,i,T)
          Pti <- ZCBond_PriceDistribution(.Object,t,alphai)          
          a <- A(alpha,beta,volDefault,alphai)          
          b <- B(beta,volDefault,alphai)          
          c <- C(eta,alphai)          
          y1 <- y1 + a*exp(b*lambdat)*c*Pti*exp(-gammat*alphai)          
        }
        
        #Calcul du second terme
        Pti <- ZCBond_PriceDistribution(.Object,t,T)
        a <- A(alpha,beta,volDefault,T)
        b <- B(beta,volDefault,T)
        c <- C(eta,T)
        y2 <- A(alpha,beta,volDefault,T)
        
        #Calcul du troisième terme
        for (i in 1:100) #100 points de calcul pour l'intégrale
        {
          delta <- (T-t)/100
          iprime <- i*delta
          Pti <- ZCBond_PriceDistribution(.Object,t,t+iprime)
          a <- A(alpha,beta,volDefault,iprime)
          b <- B(beta,volDefault,iprime)
          c <- C(eta,iprime)
          g <- G(alpha,beta,volDefault,iprime)
          h <- H(alpha,beta,volDefault,iprime)
          y3 <- y3 + delta*exp(b*lambdat)*c*Pti*(g+lambdat*h)*exp(-gammat*iprime)
        }
        
        y <- couponsRate*y1+y2+(1-omega)*y3
        
        #Ajout du coupon couru
        cc <- CouponTimeValue(couponsRate,nCoupons,t,T)
      
      return(y+cc)
      } else stop("t (starting month) must be positive and lower than maturity T")
    }
)