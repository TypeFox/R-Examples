#Models paths generation
############################################
#Scenarios Paths generation
############################################



setGeneric(
  name="customPathsGeneration",
  def = function(.Object, type=NULL)
  {
    standardGeneric("customPathsGeneration")
  }
)

setMethod(
  f="customPathsGeneration",
  signature="Scenarios",
  definition=function(.Object, type=NULL)
  { 
    horizon  <- .Object@ParamsScenarios@horizon
    nScenarios <- .Object@ParamsScenarios@nScenarios  
    
    if (!missing(type) && (type != "shortRate") && (type != "stock") && (type != "realEstate") && 
          (type !="liquiditySpread") && (type !="defaultSpread"))
    {
      stop("Parameter 'type' must be either NULL or 'shortRate', 'stock', 'realEstate', '
           liquiditySpread', 'defaultSpread'")
    }
    
    if (type == "shortRate"||type == "stock"||type == "realEstate"||missing(type))
    { 
      #Variable extraction from "Scenarios" object
      f0 <- .Object@ForwardRates
      # Vasicek Hull-White
      vol <- .Object@ParamsScenarios@vol
      k <- .Object@ParamsScenarios@k
      
      epsilon1 <- matrix(data=rnorm(nScenarios*horizon),nrow=nScenarios,ncol=horizon)
      
      shortRatePaths <- matrix(data=f0[1],nrow=nScenarios,ncol=(horizon+1))
      
        for (i in 1:horizon)
      {
        t <- (i-1)
        a <- shortRatePaths[,i]*exp(-k)+f0[i+1]-f0[i]*exp(-k)
        b <- vol^2/2*(K(t,k)^2-exp(-k)*K(t,k)^2)
        c <- sqrt(L(1,vol,k))*epsilon1[,i]
        shortRatePaths[,i+1] <- a+b+c
      }      
      .Object@shortRatePaths <- shortRatePaths
    }
      
    if (type == "stock"||missing(type))
    { 
      # Black-Scholes
      stock0 <- .Object@ParamsScenarios@stock0
      volStock <- .Object@ParamsScenarios@volStock      
      
      epsilon2 <- matrix(data=rnorm(nScenarios*horizon),nrow=nScenarios,ncol=horizon)
      
      stockPaths <- matrix(data=stock0,nrow=nScenarios,ncol=(horizon+1))
       for (i in 1:horizon)
      { # Correlation
        rho <- .Object@ParamsScenarios@rho
        stockPaths[,i+1] <- stockPaths[,i]*exp((shortRatePaths[,i+1]-volStock^2/2)+volStock*(rho*epsilon1[,i]+sqrt(1-rho^2)*epsilon2[,i]))      
      }
      .Object@stockPaths <- stockPaths
    }  
      
    if (type == "realEstate"||missing(type))
    {
      # Black-Scholes
      volRealEstate <- .Object@ParamsScenarios@volRealEstate
      realEstate0 <- .Object@ParamsScenarios@realEstate0
      
      epsilon3 <- matrix(data=rnorm(nScenarios*horizon),nrow=nScenarios,ncol=horizon)
      
      realEstatePaths <- matrix(data=realEstate0,nrow=nScenarios,ncol=(horizon+1))
      
      for (i in 1:horizon)
      {
        realEstatePaths[,i+1] <- realEstatePaths[,i]*exp((shortRatePaths[,i+1]-volRealEstate^2/2)+volRealEstate*epsilon3[,i])
        
      }
      .Object@realEstatePaths <- realEstatePaths
    }
    
    if (type == "liquiditySpread"||missing(type))
    {
      # LMN
      eta <- .Object@ParamsScenarios@eta
      liquiditySpread0 <- .Object@ParamsScenarios@liquiditySpread0
            
      epsilon4 <- matrix(data=rnorm(nScenarios*horizon),nrow=nScenarios,ncol=horizon)
      
      liquiditySpreadPaths <- matrix(data=liquiditySpread0,nrow=nScenarios,ncol=(horizon+1))
      
      for (i in 1:horizon)
      {
        liquiditySpreadPaths[,i+1] <- liquiditySpreadPaths[,i]+eta*epsilon4[,i]               
      }
      .Object@liquiditySpreadPaths <- liquiditySpreadPaths      
    }
    
    if (type == "defaultSpread"||missing(type))
    { 
      defaultSpread0 <- .Object@ParamsScenarios@defaultSpread0    
      volDefault <- .Object@ParamsScenarios@volDefault
      alpha <- .Object@ParamsScenarios@alpha
      beta <- .Object@ParamsScenarios@beta
      
      epsilon5 <- matrix(data=rnorm(nScenarios*horizon),nrow=nScenarios,ncol=horizon)
      
      #Paths matrices initialization
      defaultSpreadPaths <- matrix(data=defaultSpread0,nrow=nScenarios,ncol=(horizon+1))      
      
      for (i in 1:horizon)
      {
        a<-alpha*(beta-defaultSpreadPaths[,i])
        b<-volDefault*sqrt(defaultSpreadPaths[,i])*epsilon5[,i]
        c<-volDefault^2/4*(epsilon5[,i]^2-1)
        defaultSpreadPaths[,i+1] <- defaultSpreadPaths[,i]+a+b+c
      }
      .Object@defaultSpreadPaths <- defaultSpreadPaths
    }    
    return(.Object)
  }
  
)
