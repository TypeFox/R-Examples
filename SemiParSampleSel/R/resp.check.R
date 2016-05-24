resp.check <- function (y, margin = "N", bd = 1) {
    m2 <- c("N", "G", "P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", "LG", 
            "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", 
            "ZALG", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG")
    y <- na.omit(y)
    if (!(margin %in% m2)) 
        stop("Error in margin value. It can be: N, G, P, NB, D, PIG, S, BB, BI, GEOM, LG, 
            NBII, WARING, YULE, ZIBB, ZABB, ZABI, ZIBI, 
            ZALG, ZANBI, ZINBI, ZAP, ZIP, ZIP2, ZIPIG.")
    if (margin %in% c("G", "P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", 
            "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", 
            "ZALG", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG") && min(y) < 0) 
        stop("The variable of interest must be positive.")
    if (margin %in% c("LG") && min(y) < 1) 
        stop("The variable of interest must be greater than 1.")
    margins <- c("probit", margin)
    VC <- list(X1 = NULL, 
               X2 = matrix(1, nrow = length(y), ncol = 1), 
               dat=cbind(rep(1, length(y)), y),
               X1.d2 = 1, 
               X2.d2 = 1,
               gp1 = 1, 
               gp2 = 1,
               n = length(y),
               bd = bd,
               l.sp1 = 0, 
               l.sp2 = 0,
               l.sp3 = 0,
               l.sp4 = 0,
               l.sp5 = 0,
               gamma = 1,
               weights = 1,
               fp = FALSE,
               BivD = NULL, 
               margins = c("N", margin))
    if (margin %in% c("P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", "LG", 
                      "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", 
                      "ZALG", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG")) {
      if (margin %in% c("P")) 
        {start.val <- log(mean((y + mean(y))/2))}
      else if (margin %in% c("BI")) 
        {start.val <- qlogis(mean((y + 0.5)/(bd + 1)))}
      else if (margin %in% c("GEOM", "YULE")) 
        {start.val <- log(mean(y))}
      else if (margin %in% c("LG")) 
        {start.val <- qlogis(0.9)}
      else if (margin %in% c("NB", "PIG"))
      {start.val <- c(log(mean((y + mean(y))/2)), log(max(((var(y)/mean(y)) - 1), 0.1)))}
      else if (margin %in% c("BB"))
      {start.val <- c(qlogis(mean((y + 0.5)/(bd + 1))), log(1))}
      else if (margin %in% c("NBII"))
      {start.val <- c(log(mean((y + mean(y))/2)), log(max(((var(y)/mean(y)) - 1), 0.1)))}
      else if (margin %in% c("WARING"))
      {start.val <- c(log(mean(y)), log(1))}
      else if (margin %in% c("ZABI", "ZIBI", "ZALG"))
      {start.val <- c(qlogis(0.5), qlogis(0.3))}
      else if (margin %in% c("ZAP"))
      {start.val <- c(log(mean((y + mean(y))/2)), qlogis(0.3))}
      else if (margin %in% c("ZIP"))
      {start.val <- c(log(mean((y + mean(y))/2)), qlogis(0.1))}
      else if (margin %in% c("ZIP2"))
      {start.val <- c(log(mean((y + mean(y))/2)), qlogis(0.3))}
      else if (margin %in% c("D"))
      {start.val <- c(log(mean((y + mean(y))/2)), log(max(((var(y) - mean(y))/(mean(y)^2)), 0.1)), qlogis(0.5))}
      else if (margin %in% c("S"))
      {start.val <- c(log(mean((y + mean(y))/2)), log(max(((var(y) - mean(y))/(mean(y)^2)), 0.1)), -0.5)}
      else if (margin %in% c("ZIBB"))
      {start.val <- c(qlogis(0.5), log(0.5), qlogis(0.3))}
      else if (margin %in% c("ZABB"))
      {start.val <- c(qlogis(0.5), log(0.2), qlogis(0.3))}
      else if (margin %in% c("ZANBI"))
      {start.val <- c(log(mean((y + mean(y))/2)), log(max(((var(y) - mean(y))/(mean(y)^2)), 0.1)), qlogis(max(sum(y == 0)/length(y), 0.01)))}
      else if (margin %in% c("ZINBI"))
      {start.val <- c(log(mean((y + mean(y))/2)), log(max(((var(y) - mean(y))/(mean(y)^2)), 0.1)), qlogis(((sum(y == 0)/length(y)) + 0.01)/2))}
      else if (margin %in% c("ZIPIG"))
      {start.val <- c(log(mean((y + mean(y))/2)), log(max(((var(y) - mean(y))/(mean(y)^2)), 0.1)), qlogis(0.1))}
      ghs <- ghssDuniv
      univfit <- trust(ghs, parinit=start.val, rinit = 1, rmax = 100, dat=VC$dat, VC = VC, sp = NULL, qu.mag = NULL)
      }
    else if (margin == "N")
      {univfit <- lm(y~1)}
    else if (margin == "G")
      {univfit <- glm(y~1, family=Gamma)}
    precision <- 10^(-7)
		if (margin=="N") {
        mu <- univfit$coefficients
        sigma <- summary(univfit)$sigma
        F2 <-  pnorm(y, mean=mu, sd=sigma)
        f2 <-  dnorm(y, mean=mu, sd=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)    
        f2 <- ifelse(f2>precision, f2, precision)
     } else if (margin=="G") {
        shape <- (1/summary(univfit)$dispersion)
        rate <- mean(predict(univfit))*(1/summary(univfit)$dispersion)
     		F2 <- pgamma(y, shape=shape, rate=rate)
     		f2 <- dgamma(y, shape=shape, rate=rate) 
     		F2 <- ifelse(F2<(1-precision), F2, 1-precision)
     		F2 <- ifelse(F2>precision, F2, precision)    
     		f2 <- ifelse(f2>precision, f2, precision)
     } else  if (margin=="P") {
     		mu <- exp(univfit$argument[1])
        F2 <- pPO(y, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dPO(y, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pPO(y, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2) 
     } else if (margin=="NB") {
     		mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        F2 <- pNBI(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dNBI(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pNBI(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="D") {
      	mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        nu <- plogis(univfit$argument[3])
        F2 <- pDEL(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dDEL(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pDEL(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="PIG") {
      	mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        F2 <- pPIG(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dPIG(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pPIG(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="S") {
        mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        nu <- univfit$argument[3]
        F2 <- pSICHEL(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dSICHEL(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pSICHEL(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     }  else if (margin=="BB") {
        mu <- plogis(univfit$argument[1])
        sigma <- plogis(univfit$argument[2])
        F2 <- pBB(y, bd=bd, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dBB(y, bd=bd, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pBB(y, bd=bd, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="BI") {
        mu <- plogis(univfit$argument[1])
        F2 <- pBI(y, bd=bd, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dBI(y, bd=bd, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pBI(y, bd=bd, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="GEOM") {
        mu <- exp(univfit$argument[1])
        F2 <- pGEOM(y, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dGEOM(y, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pGEOM(y, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="LG") {
        mu <- plogis(univfit$argument[1])
        F2 <- pLG(y, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dLG(y, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pLG(y, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="NBII") {
        mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        F2 <- pNBII(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dNBII(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pNBII(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="WARING") {
        mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        nu <- univfit$argument[3]
        F2 <- pWARING(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dWARING(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pWARING(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="YULE") {
        mu <- exp(univfit$argument[1])
        F2 <- pYULE(y, mu=mu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dYULE(y, mu=mu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pYULE(y, mu=mu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZABB") {
        mu <- plogis(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        nu <- plogis(univfit$argument[3])
        F2 <- pZABB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZABB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZABB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZIBB") {
        mu <- plogis(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        nu <- plogis(univfit$argument[3])
        F2 <- pZIBB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIBB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIBB(y, bd=bd, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZABI") {
        mu <- plogis(univfit$argument[1])
        sigma <- plogis(univfit$argument[2])
        F2 <- pZABI(y, bd=bd, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZABI(y, bd=bd, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZABI(y, bd=bd, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZIBI") {
        mu <- plogis(univfit$argument[1])
        sigma <- plogis(univfit$argument[2])
        F2 <- pZIBI(y, bd=bd, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIBI(y, bd=bd, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIBI(y, bd=bd, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZALG") {
        mu <- plogis(univfit$argument[1])
        sigma <- plogis(univfit$argument[2])
        F2 <- pZALG(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZALG(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZALG(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZANBI") {
        mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        nu <- plogis(univfit$argument[3])
        F2 <- pZANBI(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZANBI(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZANBI(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZINBI") {
        mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        nu <- plogis(univfit$argument[3])
        F2 <- pZINBI(y, mu=mu, sigma=sigma, nu=nu)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZINBI(y, mu=mu, sigma=sigma, nu=nu)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZINBI(y, mu=mu, sigma=sigma, nu=nu)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZAP") {
        mu <- exp(univfit$argument[1])
        sigma <- plogis(univfit$argument[2])
        F2 <- pZAP(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZAP(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZAP(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZIP") {
        mu <- exp(univfit$argument[1])
        sigma <- plogis(univfit$argument[2])
        F2 <- pZIP(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIP(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIP(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZIP2") {
        mu <- exp(univfit$argument[1])
        sigma <- plogis(univfit$argument[2])
        F2 <- pZIP2(y, mu=mu, sigma=sigma)
        F2 <- ifelse(F2<(1-precision), F2, 1-precision)
        F2 <- ifelse(F2>precision, F2, precision)
        f2 <- dZIP2(y, mu=mu, sigma=sigma)
        f2 <- ifelse(f2>precision, f2, precision)
        F22 <- pZIP2(y, mu=mu, sigma=sigma)-f2
        F22 <- ifelse(F22<(1-precision), F22, 1-precision)
        F22 <- ifelse(F22>precision, F22, precision)
        diffs <- runif(y, F22, F2)
     } else if (margin=="ZIPIG") {
        mu <- exp(univfit$argument[1])
        sigma <- exp(univfit$argument[2])
        nu <- plogis(univfit$argument[3])
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
    par(mfrow = c(1, 2))
    hist(y, freq = FALSE, ylim = c(0, max(f2)), xlim = c(0, max(y)), main = "Histogram and Density of Response", xlab = "Response")
    xspline(sort(y), f2[order(y)], lwd = 2)
   if (margin %in% c("P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", "LG", 
                      "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", 
                      "ZALG", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG")) {
      qqnorm(qnorm(diffs), main = "Normalized Q-Q plot")
      qqline(qnorm(diffs), col = "red")
    } else if (margin == "N"){
      qqnorm(y)
      qqline(y, col = "red")
    } else if (margin == "G")  {
      qqplot(qgamma(ppoints(500), shape=shape, rate=rate), y, xlab="Theoretical Quantiles ", ylab="Sample Quantiles",  main = "Gamma Q-Q plot")
      qqline(y, distribution = function(p) qgamma(p, shape=shape, rate=rate), col = 2)
    }
  
  
 
}