#Corporate bond price distribution
#t integer
#T (year)
#Omega represents the part recovered in case of default
############################################


setGeneric(
  name="CDSPremium_PriceDistribution",
  def = function(.Object,t,T,omega)
  {
    standardGeneric("CDSPremium_PriceDistribution")
  }
)
setMethod(
  f="CDSPremium_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object,t,T,omega)
  {     
    if (!is.numeric(t)||!is.numeric(T)||!is.numeric(omega)){stop("All the parameters have to be numeric")}
    
    if (floor(t) != t){stop("Parameter t (starting month) must be an integer")}
    
    if (t == T){return (0)}
    
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
      #Initialization
      nScenarios <- .Object@ParamsScenarios@nScenarios
      y <- rep(0,nScenarios)
      y1 <- rep(0,nScenarios)
      y2 <- rep(0,nScenarios)
      gammat <- .Object@liquiditySpreadPaths[,t+1]
      lambdat <- .Object@defaultSpreadPaths[,t+1]
    }    
    
    if (t>=0 && t<=T)
    {
      #Scenarios' properties reading
      alpha <- .Object@ParamsScenarios@alpha
      beta <- .Object@ParamsScenarios@beta
      eta <- .Object@ParamsScenarios@eta
      volDefault <- .Object@ParamsScenarios@volDefault
        
      delta <- (T-t)/100
        #Numerator
        for (i in 1:100)
        { 
          iprime <- i*delta
          Pti <- ZCBond_PriceDistribution(.Object,t,t+iprime)
          b <- B(beta,volDefault,iprime)
          g <- G(alpha,beta,volDefault,iprime)
          h <- H(alpha,beta,volDefault,iprime)          
          y1 <- y1 + delta*exp(b*lambdat)*Pti*(g+lambdat*h)
        }
        
        #Denominator
        for (i in 1:100) #100 points for integral calculus
        {
          iprime <- i*delta
          Pti <- ZCBond_PriceDistribution(.Object,t,t+iprime)
          a <- A(alpha,beta,volDefault,iprime)
          b <- B(beta,volDefault,iprime)
          y2 <- y2 + delta*a*exp(b*lambdat)*Pti
        }
        
        y <- omega*y1/y2
      
      return(y)
      } else stop("t (starting month) must be positive and lower than maturity T")
    }
)