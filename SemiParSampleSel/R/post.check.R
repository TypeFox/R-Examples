post.check <- function (x, bd = 1) {
  
    sigma <- nu  <- NULL
    
    y <- x$y2[x$y1==1]
    eta2 <- x$eta2[x$y1==1]
      
  if(x$margins[2] %in% c("N", "G", "NB", "PIG", "D", "S", "BB", "NBII", "WARING", "ZIBB", "ZABB", 
                       "ZANBI", "ZINBI", "ZIPIG", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") ) { 
    sigma <- x$sigma
    if(length(x$formula) > 2) sigma <- x$sigma[x$y1==1]
  } 
  
  if(x$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") ) { nu <- x$nu
     if(length(x$formula) > 2)  nu <- x$nu[x$y1==1]   
  }
  
    precision <- 10^(-7)
		if (x$margins[2]=="N") {
        mu <- eta2
        F2 <-  pnorm(y, mean=mu, sd=sigma)
        f2 <-  dnorm(y, mean=mu, sd=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)    
        f2 <- ifelse(f2>precision, f2, precision)
     } else if (x$margins[2]=="G") {
        k <- sigma
        shape <- k
        rate <- k*exp(-eta2)
     		F2 <- pgamma(y, shape=shape, rate=rate)
     		f2 <- dgamma(y, shape=shape, rate=rate) 
     		F2 <- ifelse(F2<(1-precision), F2, 1-precision)
     		F2 <- ifelse(F2>precision, F2, precision)    
     		f2 <- ifelse(f2>precision, f2, precision)
     } else  if (x$margins[2]=="P") {
     		mu <- exp(eta2)
        F2 <- pPO(y, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dPO(y, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pPO(y, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2) 
     } else if (x$margins[2]=="NB") {
     		mu <- exp(eta2)
        F2 <- pNBI(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dNBI(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pNBI(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="D") {
      	mu <- exp(eta2)
        F2 <- pDEL(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dDEL(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pDEL(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="PIG") {
      	mu <- exp(eta2)
        F2 <- pPIG(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dPIG(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pPIG(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="S") {
        mu <- exp(eta2)
        F2 <- pSICHEL(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dSICHEL(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pSICHEL(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     }  else if (x$margins[2]=="BB") {
        mu <- plogis(eta2)
        F2 <- pBB(y, bd=bd, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dBB(y, bd=bd, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pBB(y, bd=bd, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="BI") {
        mu <- plogis(eta2)
        F2 <- pBI(y, bd=bd, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dBI(y, bd=bd, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pBI(y, bd=bd, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="GEOM") {
        mu <- exp(eta2)
        F2 <- pGEOM(y, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dGEOM(y, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pGEOM(y, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="LG") {
        mu <- plogis(eta2)
        F2 <- pLG(y, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dLG(y, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pLG(y, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="NBII") {
        mu <- exp(eta2)
        F2 <- pNBII(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dNBII(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pNBII(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="WARING") {
        mu <- exp(eta2)
        F2 <- pWARING(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dWARING(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pWARING(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="YULE") {
        mu <- exp(eta2)
        F2 <- pYULE(y, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dYULE(y, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pYULE(y, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZABB") {
        mu <- plogis(eta2)
        F2 <- pZABB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZABB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZABB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZIBB") {
        mu <- plogis(eta2)
        F2 <- pZIBB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIBB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIBB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZABI") {
        mu <- plogis(eta2)
        F2 <- pZABI(y, bd=bd, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZABI(y, bd=bd, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZABI(y, bd=bd, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZIBI") {
        mu <- plogis(eta2)
        F2 <- pZIBI(y, bd=bd, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIBI(y, bd=bd, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIBI(y, bd=bd, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZALG") {
        mu <- plogis(eta2)
        F2 <- pZALG(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZALG(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZALG(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZANBI") {
        mu <- exp(eta2)
        F2 <- pZANBI(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZANBI(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZANBI(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZINBI") {
        mu <- exp(eta2)
        F2 <- pZINBI(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZINBI(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZINBI(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZAP") {
        mu <- exp(eta2)
        F2 <- pZAP(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZAP(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZAP(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZIP") {
        mu <- exp(eta2)
        F2 <- pZIP(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIP(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIP(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZIP2") {
        mu <- exp(eta2)
        F2 <- pZIP2(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIP2(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIP2(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (x$margins[2]=="ZIPIG") {
        mu <- exp(eta2)
        F2 <- pZIPIG(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIPIG(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIPIG(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     }               
    #par(mfrow = c(1, 1))
    #par(mfrow = c(1, 2))
    #hist(y, freq = FALSE, ylim = c(0, max(f2)), xlim = c(0, max(y)), breaks = breaks, main = "Histogram and Density of Response", xlab = "Response")
    #  if(kernel.dens == TRUE) {  xspline(sort(y), f2[order(y)], lwd = 2)  }
   if (x$margins[2] %in% c("P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", "LG", 
                      "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", 
                      "ZALG", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG")) {
      qqnorm(qnorm(diffs), main = "Normalized Q-Q plot")
      qqline(qnorm(diffs), col = "red")
    } else if (x$margins[2] == "N"){
      qqnorm(y)
      qqline(y, col = "red")
    } else if (x$margins[2] == "G")  {
      qqplot(qgamma(ppoints(500), shape=shape, rate=rate), y, xlab="Theoretical Quantiles ", ylab="Sample Quantiles",  main = "Gamma Q-Q plot")
      qqline(y, distribution = function(p) qgamma(p, shape=shape, rate=rate), col = 2)
    }
  
  
 
}
