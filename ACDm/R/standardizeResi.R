standardizeResi <- function(fitModel, transformation = "probIntegral"){ 
  
  if(!("acdFit" %in% class(fitModel))) stop("fitModel is not of class 'acdFit'")
  
  transformation <- match.arg(transformation, c("probIntegral", "cox-snell"))
  
  if(fitModel$distribution == "exponential"){ 
    
    returnValue <- stats::pexp(fitModel$residuals)
    
  } else if(fitModel$distribution == "weibull"){            
    
    returnValue <- stats::pweibull(fitModel$residuals, shape = fitModel$dPara, scale = fitModel$forcedDistPara)
    
  } else if(fitModel$distribution == "burr"){
    
    returnValue <- pburr(fitModel$residuals, theta = fitModel$forcedDistPara, kappa = fitModel$dPara[1], sig2 = fitModel$dPara[2])
    
  } else if(fitModel$distribution == "gengamma"){
    
    returnValue <- pgengamma(fitModel$residuals, gamma = fitModel$dPara[2], kappa = fitModel$dPara[1], lambda = fitModel$forcedDistPara)
    
  } else if(fitModel$distribution == "qweibull"){
    
    returnValue <- pqweibull(fitModel$residuals, a = fitModel$dPara[1], qdist = fitModel$dPara[2], b = fitModel$forcedDistPara)
    
  } else if(fitModel$distribution == "mixqwe"){
    
    returnValue <- pmixqwe(fitModel$residuals, pdist = fitModel$dPara[1], a = fitModel$dPara[2], qdist = fitModel$dPara[3], lambda = fitModel$dPara[4], b = fitModel$forcedDistPara)
    
  } else if(fitModel$distribution == "mixqww"){
    
    returnValue <- pmixqww(fitModel$residuals, pdist = fitModel$dPara[1], a = fitModel$dPara[2], qdist = fitModel$dPara[3], theta = fitModel$dPara[4], gamma = fitModel$dPara[5], b = fitModel$forcedDistPara)
    
  } else if(fitModel$distribution == "mixinvgauss"){
    
    returnValue <- pmixinvgauss(fitModel$residuals, theta = fitModel$dPara[1], lambda = fitModel$dPara[2], gamma = fitModel$dPara[3])
    
  } else stop("not yet implemented for the ", fitModel$distribution, " distribution")
  
  if(transformation == "cox-snell"){
    if(fitModel$distribution == "exponential"){ 
      
      warning("Cox-Snell transformation for the exponential distribution will leave the residuals unchanged!")
      return(fitModel$residuals)
      
    } else{            
      
      return(-log(1 - returnValue))
      
    } 
  } else returnValue
  
}