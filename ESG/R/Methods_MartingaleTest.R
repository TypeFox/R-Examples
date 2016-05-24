#Martingale test
############################################

setGeneric(
  name="MartingaleTest",
  def = function(.Object)
  {
    standardGeneric("MartingaleTest")
  }
)
setMethod(
  f="MartingaleTest",
  signature="Scenarios",
  definition=function(.Object)
  {
    nScenarios <- .Object@ParamsScenarios@nScenarios
    horizon <- .Object@ParamsScenarios@horizon
    stock0 <- .Object@ParamsScenarios@stock0
    
    #Calculus support matrices
    tool1 <- matrix(data=1,nrow=(horizon+1),ncol=(horizon+1))
    tool1[1,] <- 0
    tool1[lower.tri(tool1)] <- 0
    tool2 <- rep(1/nScenarios,nScenarios)
    
    
    rateSum <- .Object@shortRatePaths%*%tool1
    actuFactors <- exp(-rateSum)
    result <- .Object@stockPaths*actuFactors
    
    mean <- tool2%*%result/stock0
    
    plot(0:horizon,mean)
    lines(0:horizon,rep(1,horizon+1))
    title(main = "Annual mean of discounted prices")
    return(max(abs(1-mean)))
  }
)