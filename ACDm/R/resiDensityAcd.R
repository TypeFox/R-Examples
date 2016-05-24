resiDensityAcd <- function(fitModel, xlim = NULL, binwidth = .1, density = FALSE){
  
  
  ..density.. <- dexp <- dweibull <- residuals <- NULL
  
  df <- data.frame(residuals = fitModel$residuals)
  
  #sets the PDF function for the assumed distribution:
  distFnc <- switch(fitModel$distribution,
                    "exponential" = dexp,
                    "weibull" = dweibull,
                    "burr" = dburr,
                    "gengamma" = dgengamma,
                    "genf" = dgenf,
                    "qweibull" = dqweibull,
                    "mixqwe" = dmixqwe,
                    "mixqww" = dmixqww,
                    "mixinvgauss" = dmixinvgauss)
  
  if(fitModel$distribution == "exponential"){
    
    paraList <- list(rate = 1)
    
  } else if(fitModel$distribution == "weibull"){ 
    
    paraList <- list(shape = fitModel$dPara, scale = 1/(gamma(1+1/fitModel$dPara)))
    
  } else if(fitModel$distribution == "burr"){ 
    
    kappa = fitModel$dPara[1]
    sig2 = fitModel$dPara[2]
    
    theta <- ((gamma(1+1/kappa)*gamma(1/sig2 - 1/kappa))/(sig2^(1+1/kappa)*gamma(1/sig2+1)))^(kappa)
    
    paraList <- list(kappa = kappa, sig2 = sig2, theta = theta)
      
  } else if(fitModel$distribution == "gengamma"){ 
    
    kappa <- fitModel$dPara[1]
    gamma <- fitModel$dPara[2]
    lambda <- fitModel$forcedDistPara
    
    paraList <- list(gamma = gamma, kappa = kappa, lambda = lambda)
    
  } else if(fitModel$distribution == "genf"){ 
    
    kappa <- fitModel$dPara[1]
    eta <- fitModel$dPara[2]
    gamma <- fitModel$dPara[3]
    lambda <- fitModel$forcedDistPara
    
    paraList <- list(kappa = kappa, eta = eta, gamma = gamma, lambda = lambda)
    
  } else if(fitModel$distribution == "qweibull"){ 
    
    a <- fitModel$dPara[1]
    q <- fitModel$dPara[2]
    b <- fitModel$forcedDistPara      
    
    paraList <- list(a = a, q = q, b = b)
    
  } else if(fitModel$distribution == "mixqwe"){ 
    
    p <- fitModel$dPara[1]
    a <- fitModel$dPara[2]
    qdist <- fitModel$dPara[3]
    lambda <- fitModel$dPara[4]
    b <- fitModel$forcedDistPara    
    
    paraList <- list(pdist = p, a = a, qdist = qdist, lambda = lambda, b = b)
    
  } else if(fitModel$distribution == "mixqww"){ 
    
    p <- fitModel$dPara[1]
    a <- fitModel$dPara[2]
    qdist <- fitModel$dPara[3]
    theta <- fitModel$dPara[4]
    gamma <- fitModel$dPara[5]
    b <- fitModel$forcedDistPara
    
    paraList <- list(pdist = p, a = a, qdist = qdist, theta = theta, gamma = gamma, b = b)
    
  } else if(fitModel$distribution == "mixinvgauss"){ 
    
    theta <- fitModel$dPara[1]
    lambda <- fitModel$dPara[2]
    gamma <- fitModel$dPara[3]
    
    paraList <- list(theta = theta, lambda = lambda, gamma = gamma, forceExpectation = T)
    
  } else stop("the provided distribution does not exist")  
  
  g <- ggplot(df, aes(x=residuals)) 
  g <- g + geom_histogram(aes(y = ..density..), binwidth = binwidth, alpha = 0.4)   
  if(density) g <- g + stat_density(aes(colour = 'Empirical'), cex = 1, geom = "line", adjust = .2)    
  g <- g + stat_function(fun = distFnc, aes(colour = 'Implied'), 
                         args = paraList, n = 5000, cex = 1)
  if(length(xlim) != 0) g <- g + xlim(xlim)
  g + scale_colour_manual(name = 'Density', values = c('red', 'blue'))
  
}