##
## Augmented-Dickey-Fuller Test
##
ur.df <- function (y, type = c("none", "drift", "trend"), lags = 1, selectlags = c("Fixed", "AIC", "BIC")) 
{
    selectlags<-match.arg(selectlags)
    type <- match.arg(type)
    if (ncol(as.matrix(y)) > 1) 
        stop("\ny is not a vector or univariate time series.\n")
    if (any(is.na(y))) 
        stop("\nNAs in y.\n")
    y <- as.vector(y)
    lag <- as.integer(lags)
    if (lag < 0) 
        stop("\nLags must be set to an non negative integer value.\n")
    CALL <- match.call()
    DNAME <- deparse(substitute(y))
    x.name <- deparse(substitute(y))
    lags <- lags + 1
    z <- diff(y)
    n <- length(z)
    x <- embed(z, lags)
    z.diff <- x[, 1]
    z.lag.1 <- y[lags:n]
    tt <- lags:n
    if (lags > 1) {
    if(selectlags!="Fixed"){
      critRes<-rep(NA, lags)
      for(i in 2:(lags)){
          z.diff.lag = x[, 2:i]
	  if (type == "none") 
            result <- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
	  if (type == "drift") 
            result <- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)  
	  if (type == "trend") 
            result <- lm(z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
	  critRes[i]<-AIC(result, k = switch(selectlags, "AIC" = 2, "BIC" = log(length(z.diff))))
	}
	lags<-which.min(critRes)
    }
        z.diff.lag = x[, 2:lags]
        if (type == "none") {
            result <- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
            tau <- coef(summary(result))[1, 3]
            teststat <- as.matrix(tau)
            colnames(teststat) <- 'tau1'
          }
        if (type == "drift") {
            result <- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
            tau <- coef(summary(result))[2, 3]
            phi1.reg <- lm(z.diff ~ -1 + z.diff.lag)
            phi1 <- anova(phi1.reg, result)$F[2]
            teststat <- as.matrix(t(c(tau, phi1)))
            colnames(teststat) <- c('tau2', 'phi1')
          }
        if (type == "trend") {
            result <- lm(z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
            tau <- coef(summary(result))[2, 3]
            phi2.reg <- lm(z.diff ~ -1 + z.diff.lag)
            phi3.reg <- lm(z.diff ~ z.diff.lag)
            phi2 <- anova(phi2.reg, result)$F[2]
            phi3 <- anova(phi3.reg, result)$F[2]
            teststat <- as.matrix(t(c(tau, phi2, phi3)))
            colnames(teststat) <- c('tau3', 'phi2', 'phi3')
          }
    }
    else {
        if (type == "none") {
            result <- lm(z.diff ~ z.lag.1 - 1)
            tau <- coef(summary(result))[1, 3]
            teststat <- as.matrix(tau)
            colnames(teststat) <- 'tau1'
        }
        if (type == "drift") {
            result <- lm(z.diff ~ z.lag.1 + 1)
            phi1.reg <- lm(z.diff ~ -1)
            phi1 <- anova(phi1.reg, result)$F[2]
            tau <- coef(summary(result))[2, 3]
            teststat <- as.matrix(t(c(tau, phi1)))
            colnames(teststat) <- c('tau2', 'phi1')
        }
        if (type == "trend") {
            result <- lm(z.diff ~ z.lag.1 + 1 + tt)
            phi2.reg <- lm(z.diff ~ -1)
            phi3.reg <- lm(z.diff ~ 1)
            phi2 <- anova(phi2.reg, result)$F[2]
            phi3 <- anova(phi3.reg, result)$F[2]
            tau <- coef(summary(result))[2, 3]
            teststat <- as.matrix(t(c(tau, phi2, phi3)))
            colnames(teststat) <- c('tau3', 'phi2', 'phi3')
        }
    }
    rownames(teststat) <- 'statistic'
    testreg <- summary(result)
    res <- residuals(testreg)
    if(n < 25)
      rowselec <- 1
    if(25 <= n & n < 50)
      rowselec <- 2
    if(50 <= n & n < 100)
      rowselec <- 3
    if(100 <= n & n < 250)
      rowselec <- 4
    if(250 <= n & n < 500)
      rowselec <- 5
    if(n >= 500)
      rowselec <- 6
    if (type == "none"){ 
        cval.tau1 <- rbind(
                           c(-2.66, -1.95, -1.60),
                           c(-2.62, -1.95, -1.61),
                           c(-2.60, -1.95, -1.61),
                           c(-2.58, -1.95, -1.62),
                           c(-2.58, -1.95, -1.62),
                           c(-2.58, -1.95, -1.62))
        cvals <- t(cval.tau1[rowselec, ])
        testnames <- 'tau1'
      }
    if (type == "drift"){ 
        cval.tau2 <- rbind(
                           c(-3.75, -3.00, -2.63),
                           c(-3.58, -2.93, -2.60),
                           c(-3.51, -2.89, -2.58),
                           c(-3.46, -2.88, -2.57),
                           c(-3.44, -2.87, -2.57),
                           c(-3.43, -2.86, -2.57))
        cval.phi1 <- rbind(
                           c(7.88, 5.18, 4.12),
                           c(7.06, 4.86, 3.94),
                           c(6.70, 4.71, 3.86),
                           c(6.52, 4.63, 3.81),
                           c(6.47, 4.61, 3.79),
                           c(6.43, 4.59, 3.78))
        cvals <- rbind(
                      cval.tau2[rowselec, ],
                      cval.phi1[rowselec, ])
        testnames <- c('tau2', 'phi1')
      }
    if (type == "trend"){ 
        cval.tau3 <- rbind(
                           c(-4.38, -3.60, -3.24),
                           c(-4.15, -3.50, -3.18),
                           c(-4.04, -3.45, -3.15),
                           c(-3.99, -3.43, -3.13),
                           c(-3.98, -3.42, -3.13),
                           c(-3.96, -3.41, -3.12))
        cval.phi2 <- rbind(
                           c(8.21, 5.68, 4.67),
                           c(7.02, 5.13, 4.31),
                           c(6.50, 4.88, 4.16),
                           c(6.22, 4.75, 4.07),
                           c(6.15, 4.71, 4.05),
                           c(6.09, 4.68, 4.03))
        cval.phi3 <- rbind(
                           c(10.61, 7.24, 5.91),
                           c( 9.31, 6.73, 5.61),
                           c( 8.73, 6.49, 5.47),
                           c( 8.43, 6.49, 5.47),
                           c( 8.34, 6.30, 5.36),
                           c( 8.27, 6.25, 5.34))  
        cvals <- rbind(
                      cval.tau3[rowselec, ],
                      cval.phi2[rowselec, ],
                      cval.phi3[rowselec, ])

        testnames <- c('tau3', 'phi2', 'phi3')
      }
    colnames(cvals) <- c("1pct", "5pct", "10pct")
    rownames(cvals) <- testnames
   
    new("ur.df", y = y, model = type, cval=cvals, lags=lag, teststat = teststat, testreg=testreg, res=res, test.name="Augmented Dickey-Fuller Test")
}
