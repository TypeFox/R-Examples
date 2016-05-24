plotHazard <- function(fitModel, breaks = 20, implied = TRUE, xstop = NULL){
    
  hazard <- residual <- errorTerm <- NULL
    
  if("acdFit" %in% class(fitModel)){
    e <- stats::quantile(fitModel$residuals, seq(0, 1 - 1/breaks, 1/breaks))
    h <- (1/breaks)/(1-1:(breaks-1)/breaks+1/(2*breaks))*(1/(e[2:breaks]-e[1:(breaks-1)]))
    e2 <- stats::quantile(fitModel$residuals, seq(0 + 1/(2*breaks), 1 - 3/(2*breaks), 1/breaks))
    df <- data.frame(residual = e2, hazard = h, curve = "Nonparametric")
    
    g <- ggplot(df, aes(y=hazard, x=residual))
    g <- g + geom_line() + geom_point() + ylab("hazard") + xlab("residual")
    
    if(length(xstop) == 0) xstop <- max(e2)+1
    xstart <- 0.01
    
    if(fitModel$distribution == "weibull"){
      
      gamma <- fitModel$dPara
      theta <- .returnFixedMeanPara(2, gamma)
            
      df2 <- data.frame(errorTerm = seq(xstart, xstop,.01), hazard = theta * gamma * seq(xstart, xstop,.01)^(gamma-1), curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2)  + ggtitle("Hazard function estimates: nonparametric (black) and Weibull implied (red).") 
      
      print(g)
      
    } else if(fitModel$distribution == "burr"){
      
      kappa <- fitModel$dPara[1]
      sig2 <- fitModel$dPara[2]   
      theta <- ((gamma(1+1/kappa)*gamma(1/sig2 - 1/kappa))/(sig2^(1+1/kappa)*gamma(1/sig2+1)))^(kappa)
      
      df2 <- data.frame(errorTerm = seq(xstart, xstop,.01), hazard = (theta*kappa*seq(xstart, xstop,.01)^(kappa-1))/(1+sig2*theta*seq(xstart, xstop,.01)^(kappa)), curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2) + ggtitle("Hazard function estimates: nonparametric (black) and Burr implied (red).") 
      
      print(g)
      
    } else if(fitModel$distribution == "exponential"){
      
      df2 <- data.frame(errorTerm = seq(xstart, xstop,.1), hazard = 1, curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2) + ggtitle("Hazard function estimates: nonparametric (black) and exponential implied (red).") 
      
      print(g)
      
    } else if(fitModel$distribution == "gengamma"){
      
      kappa <- fitModel$dPara[1]
      gammaPara <- fitModel$dPara[2]
      df2 <- data.frame(errorTerm = seq(xstart, xstop,.01), 
                        hazard = gengammaHazard(seq(xstart, xstop,.01), gamma = gammaPara, kappa = kappa, forceExpectation = T), curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2) + ggtitle("Hazard function estimates: nonparametric (black) and generelized Gamma implied (red).") 
      
      print(g)
      
    } else if(fitModel$distribution == "genf"){      
      
      kappa <- fitModel$dPara[1]
      eta <- fitModel$dPara[2]
      gamma <- fitModel$dPara[3]
      lambda <- fitModel$forcedDistPara      
      
      df2 <- data.frame(errorTerm = seq(xstart, xstop,.01), 
                        hazard = genfHazard(seq(xstart, xstop,.01), kappa = kappa, eta = eta, gamma = gamma, lambda = lambda), curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2) + ggtitle("Hazard function estimates: nonparametric (black) and generelized F implied (red).") 
      
      print(g)
      
    } else if(fitModel$distribution == "qweibull"){  
      
      a <- fitModel$dPara[1]
      q <- fitModel$dPara[2]
      b <- fitModel$forcedDistPara      
      
      
      df2 <- data.frame(errorTerm = seq(xstart, xstop, .01), 
                        hazard = qweibullHazard(seq(xstart, xstop,.01), a = a, qdist = qdist, b = b), curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2) + ggtitle("Hazard function estimates: nonparametric (black) and q-Weibull implied (red).") 
      
      print(g)
      
    } else if(fitModel$distribution == "mixqwe"){  
      
      p <- fitModel$dPara[1]
      a <- fitModel$dPara[2]
      qdist <- fitModel$dPara[3]
      lambda <- fitModel$dPara[4]
      b <- fitModel$forcedDistPara
      
      df2 <- data.frame(errorTerm = seq(xstart, xstop, .01), 
                        hazard = mixqweHazard(seq(xstart, xstop,.01), pdist = p, a = a, qdist = qdist, lambda = lambda, b = b), curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2) + ggtitle("Hazard function estimates: nonparametric (black) and q-Weibull mixed with exponential implied (red).") 
      
      print(g)
    } else if(fitModel$distribution == "mixqww"){  
      
      p <- fitModel$dPara[1]
      a <- fitModel$dPara[2]
      qdist <- fitModel$dPara[3]
      theta <- fitModel$dPara[4]
      gamma <- fitModel$dPara[5]
      b <- fitModel$forcedDistPara
      
      df2 <- data.frame(errorTerm = seq(xstart, xstop, .01), 
                        hazard = mixqwwHazard(seq(xstart, xstop,.01), pdist = p, a = a, qdist = qdist, theta = theta, gamma = gamma, b = b), curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2) + ggtitle("Hazard function estimates: nonparametric (black) and q-Weibull mixed with Weibull implied (red).") 
      
      print(g)
      
    } else if(fitModel$distribution == "mixinvgauss"){  
      
      theta <- fitModel$dPara[1]
      lambda <- fitModel$dPara[2]
      gamma <- fitModel$dPara[3]
      
      df2 <- data.frame(errorTerm = seq(xstart, xstop, .01), 
                        hazard = mixinvgaussHazard(seq(xstart, xstop,.01), theta = theta, lambda = lambda, gamma = gamma, forceExpectation = T), curve = "Implied")
      g <- g + geom_line(data = df2, aes(y=hazard, x=errorTerm), linetype = 1, colour = 2) + ggtitle("Hazard function estimates: nonparametric (black) and finite mixature of inverse Gaussian implied (red).") 
      
      print(g)
      
    } else{
      g <- g + ggtitle("Hazard function estimates: nonparametric") 
      print(g)
    }
  } else if (is.numeric(fitModel)){
    e <- stats::quantile(fitModel, seq(0, 1 - 1/breaks, 1/breaks))
    h <- (1/breaks)/(1-1:(breaks-1)/breaks+1/(2*breaks))*(1/(e[2:breaks]-e[1:(breaks-1)]))
    e2 <- stats::quantile(fitModel, seq(0 + 1/(2*breaks), 1 - 3/(2*breaks), 1/breaks))
    df <- data.frame(residual = e2, hazard = h, curve = "Nonparametric")
    
    g <- ggplot(df, aes(y=hazard, x=residual))
    g <- g + geom_line() + geom_point() + ylab("hazard") + xlab("residual")
    g <- g + ggtitle("Hazard function estimate") 
    print(g)
  } else warning("The 'fitModel' argument has to be either an estimated model of class 'acdFit' or a numeric vector")
}