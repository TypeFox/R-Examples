# $Id: ci.R 2060 2015-07-19 03:22:30Z warnes $

ci  <-  function(x, confidence=0.95,alpha=1-confidence,...)
  UseMethod("ci")

ci.numeric <- function(x, confidence=0.95,alpha=1-confidence,na.rm=FALSE,...)
  {
    est <- mean(x, na.rm=na.rm)
    stderr <-  sd(x, na.rm=na.rm)/sqrt(nobs(x));
    ci.low  <- est + qt(alpha/2, nobs(x)-1) * stderr
    ci.high <- est - qt(alpha/2, nobs(x)-1) * stderr
    retval  <- c(
                 Estimate=est,
                 "CI lower"=ci.low,
                 "CI upper"=ci.high,
                 "Std. Error"=stderr
                 )

    retval
  }

ci.binom <- function(x, confidence=0.95,alpha=1-confidence,...)
  {
      if( !(all(x %in% c(0,1))) ) stop("Binomial values must be either 0 or 1.")
      if( all(x==0) || all(x==1) )
          warning("All observed values are ", as.numeric(x[1]), ", so estimated Std. Error is 0.")

    est  <-  mean(x, na.rm=TRUE)
    n <- nobs(x)
    x <- sum(x)
    stderr <- sqrt(est*(1-est)/n)

    ci.low  <- qbeta(   alpha/2, x  , n + 1 - x)
    ci.high <- qbeta(1- alpha/2, x+1, n-x      )

    retval  <- cbind(Estimate=est,
                     "CI lower"=ci.low,
                     "CI upper"=ci.high,
                     "Std. Error"= stderr
                     )
    retval
  }

ci.lm  <-  function(x,confidence=0.95,alpha=1-confidence,...)
  {
    x  <-  summary(x)
    est  <-  coef(x)[,1] ;
    ci.low  <- est + qt(alpha/2, x$df[2]) * coef(x)[,2] ;
    ci.high <- est - qt(alpha/2, x$df[2]) * coef(x)[,2] ;
    retval  <- cbind(Estimate=est,
                     "CI lower"=ci.low,
                     "CI upper"=ci.high,
                     "Std. Error"= coef(x)[,2],
                     "p-value" = coef(x)[,4])

    retval
  }

ci.lme <- function(x,confidence=0.95,alpha=1-confidence,...)
  {
    x  <-  summary(x)
    est  <-  x$tTable[,"Value"] ;
    ci.low  <- est + qt(alpha/2, x$tTable[,"DF"]) * x$tTable[,"Std.Error"] ;
    ci.high <- est - qt(alpha/2, x$tTable[,"DF"]) * x$tTable[,"Std.Error"] ;
    retval  <- cbind(Estimate=est,
                     "CI lower"=ci.low,
                     "CI upper"=ci.high,
                     "Std. Error"= x$tTable[,"Std.Error"],
                     "DF" = x$tTable[,"DF"],
                     "p-value" = x$tTable[,"p-value"])
    rownames(retval)  <-  rownames(x$tTable)
    retval
  }

## ci.mer <- function (x,
##                     confidence = 0.95,
##                     alpha = 1 - confidence,
##                     n.sim = 1e4,
##                     ...)
## {
##     x.effects <- x@fixef
##     n <- length(x.effects)

##     retval <- gmodels::est.mer(obj = x,
##                                 cm = diag(n),
##                                 beta0 = rep(0, n),
##                                 conf.int = confidence,
##                                 show.beta0 = FALSE,
##                                 n.sim = n.sim)

##     retval <- retval[,
##                      c("Estimate", "Lower.CI", "Upper.CI", "Std. Error", "p value"),
##                      drop=FALSE
##                      ]
##     colnames(retval)[c(2:3, 5)] <- c("CI lower", "CI upper", "p-value")
##     rownames(retval) <- names(x.effects)

##     retval
## }


ci.estimable  <-  function(x,confidence=0.95,alpha=1-confidence,...)
  {
    ci.low  <- x$Estimate + qt(alpha/2, x$DF) * x$"Std. Error"
    ci.high <- x$Estimate - qt(alpha/2, x$DF) * x$"Std. Error"
    retval  <- cbind(Estimate=x$Estimate,
                     "CI lower"=ci.low,
                     "CI upper"=ci.high,
                     "Std. Error"= x$"Std. Error",
                     "p-value" = x$"Pr(>|t|)"
                     )
    rownames(retval) <- rownames(x)

    retval
  }
