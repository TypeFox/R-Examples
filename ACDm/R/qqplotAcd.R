qqplotAcd <- function(fitModel, xlim = NULL, ylim = NULL){
  residuals <- NULL
  
  df <- data.frame(residuals = fitModel$residuals)
  if(fitModel$distribution == "exponential"){
    
    g <- ggplot(df, aes(sample=residuals)) + stat_qq(distribution = stats::qexp, geom="point")
    if(length(xlim) != 0) g <- g + xlim(xlim)
    if(length(ylim) != 0 ) g <- g + ylim(ylim) 
    g + geom_abline(intercept = 0, slope = 1, color="red") + xlab(paste(fitModel$distribution, "theoretical quantiles"))
    
  } else if(fitModel$distribution == "weibull"){
    
    g <- ggplot(df, aes(sample=residuals)) 
    g <- g + stat_qq(distribution = stats::qweibull, dparams = list(shape = fitModel$dPara, scale = 1/(gamma(1+1/fitModel$dPara))))
    g <- g + geom_abline(intercept = 0, slope = 1, color="red") + xlab(paste(fitModel$distribution, "theoretical quantiles"))  
    if(length(xlim) != 0) g <- g + xlim(xlim)
    if(length(ylim) != 0 ) g <- g + ylim(ylim) 
    g
    
  } else if(fitModel$distribution == "burr"){
    
    burrQ <- function(p, kappa, sig2){
      theta <- ((gamma(1+1/kappa)*gamma(1/sig2 - 1/kappa))/(sig2^(1+1/kappa)*gamma(1/sig2+1)))^(kappa)
      return((((1-p)^(-sig2)-1)/(sig2*theta))^(1/kappa))
    } 
    
    g <- ggplot(df, aes(sample=residuals)) + stat_qq(distribution = burrQ, dparams = list(kappa = fitModel$dPara[1], sig2 = fitModel$dPara[2])) + geom_abline(intercept = 0, slope = 1, color="red") + xlab(paste(fitModel$distribution, "theoretical quantiles"))
    if(length(xlim) != 0) g <- g + xlim(xlim)
    if(length(ylim) != 0 ) g <- g + ylim(ylim) 
    g
    
  } else if(fitModel$distribution == "gengamma"){
    kappa <- fitModel$dPara[1]
    gammaPara <- fitModel$dPara[2]
    
    g <- ggplot(df, aes(sample=residuals)) + stat_qq(distribution = qgengamma, dparams = list(gamma = gammaPara, kappa = kappa, forceExpectation = T)) + geom_abline(intercept = 0, slope = 1, color="red") + xlab(paste(fitModel$distribution, "theoretical quantiles"))
    if(length(xlim) != 0) g <- g + xlim(xlim)
    if(length(ylim) != 0 ) g <- g + ylim(ylim) 
    g
    
  } else if(fitModel$distribution == "qweibull"){
    
    a <- fitModel$dPara[1]
    qdist <- fitModel$dPara[2]
    b <- fitModel$forcedDistPara   
    
    g <- ggplot(df, aes(sample=residuals)) + stat_qq(distribution = qqweibull, dparams = list(a = a, qdist = qdist, b = b)) + geom_abline(intercept = 0, slope = 1, color="red") + xlab(paste(fitModel$distribution, "theoretical quantiles"))
    if(length(xlim) != 0) g <- g + xlim(xlim)
    if(length(ylim) != 0 ) g <- g + ylim(ylim) 
    g
    
  } else stop("The QQ plot function is not yet implemented for this distribution")
}