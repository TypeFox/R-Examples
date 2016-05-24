######
## Weighted.Box.test -- The Weighted Box-Pierce, Ljung-Box, Monti, McLeod-Li type test for fitted ARMA and detection of nonlinear Models
##
## Quick intro for the input parameters
##
##      x       -- the residuals or initial unfitted time series
##     lag      -- The lag we wish to test at, this is m or M in the literature, default=1
##     type     -- The type of test, default="Box-Pierce"
##    fitdf     -- The number of ARMA parameter that have been fit to the series x, default=0
##
## The following are boolean flags to perform transformations to detect nonlinear processes, all default=FALSE
##
##    sqrt.res  -- Should the residuals be squared to detect for nonlinear models?  If so, fitdf is ignore, the type "Ljung-Box" is really the "McLeod-Li" type test
## log.sqrd.res -- Take the log of the squared residuals, generally simulations show this is less powerful than the sqrd.res or abs.res
##    abs.res   -- Take the absolute value of the residuals, similar to the above two boolean flags.  Simulations indicate abs.res is more powerful at detecting long memory processes
##
## For backward compatability
##
##   weighted   -- For backward compatability, If set to FALSE, you will perform the traditional Box-Pierce, Ljung-Box or Monti or McLeod Li type if using a transformation, default=TRUE.
######

Weighted.Box.test <- function (x, lag = 1, type = c("Box-Pierce", "Ljung-Box", "Monti"), fitdf = 0, sqrd.res = FALSE, log.sqrd.res = FALSE, abs.res = FALSE, weighted=TRUE)
{
 ### Error Checking
 ###
    if (NCOL(x) > 1) 
       stop("x is not a vector or univariate time series");
    if (lag < 1)
       stop("Lag must be positive");
    if (fitdf < 0)
       stop("Fitdf cannot be negative");
    if (fitdf >= lag)
       stop("Lag must exceed fitted degrees of freedom");
    if( (sqrd.res && log.sqrd.res) || (sqrd.res && abs.res) || (log.sqrd.res && abs.res) )
       stop("Only one option of: sqrd.res, log.sqrd.res or abs.res can be selected");

 ### Find the type
    DNAME <- deparse(substitute(x))
    type <- match.arg(type)

 ### Transform if checking for nonlinear
    if(abs.res) {
       x <- abs(x);
    }
    if(sqrd.res || log.sqrd.res) {
       x <- x^2;
    }
    if(log.sqrd.res) {
       x <- log(x);
    }

 if(weighted) 
 {
    if( type == "Monti") {
       METHOD <- "Weighted Monti test (Gamma Approximation)"
       cor <- acf(x, lag.max = lag, type="partial", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[1:lag];
    }
    else {
       cor <- acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[2:(lag + 1)];
    }

    if(type == "Ljung-Box") {
       METHOD <- "Weighted Ljung-Box test (Gamma Approximation)"
    }

    n <- sum(!is.na(x))

    weights <- (lag - 1:lag+1)/(lag);

    if (type == "Box-Pierce") {
       METHOD <- "Weighted Box-Pierce test (Gamma Approximation)"
       STATISTIC <- n * sum(weights*obs^2)
    }
    else {
       STATISTIC <- n * (n + 2) * sum(weights*(1/seq.int(n - 1, n - lag)*obs^2))
    }

    if(sqrd.res) {
       fitdf <- 0;
       names(STATISTIC) <- "Weighted X-squared on Squared Residuals for detecting nonlinear processes"
    }
    else if(log.sqrd.res) { 
       fitdf <- 0;
       names(STATISTIC) <- "Weighted X-squared on Log-Squared Residuals for detecting nonlinear processes";
    }
    else if(abs.res) { 
       fitdf <- 0;
       names(STATISTIC) <- "Weighted X-squared on Absolute valued Residuals for detecting nonlinear processes";
    }
    else {
       names(STATISTIC) <- "Weighted X-squared on Residuals for fitted ARMA process"
    }

    shape <- (3/4)*(lag+1)^2*lag/(2*lag^2 + 3*lag + 1 - 6*lag*fitdf);
    scale <- (2/3)*(2*lag^2 + 3*lag + 1 - 6*lag*fitdf)/lag/(lag+1);

    PARAMETER <- c(shape, scale);
    names(PARAMETER) <- c("Shape", "Scale")

    PVAL <- 1 - pgamma(STATISTIC, shape=shape, scale=scale)
    names(PVAL) <- "Approximate p-value"
 }
 else #Not weighted
 {
    if( type == "Monti") {
       METHOD <- "Monti test"
       cor <- acf(x, lag.max = lag, type="partial", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[1:lag];
    }
    else {
       cor <- acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[2:(lag + 1)];
    }

    if(type == "Ljung-Box") {
       METHOD <- "Ljung-Box test"
    }

    n <- sum(!is.na(x))


    if (type == "Box-Pierce") {
       METHOD <- "Box-Pierce test"
       STATISTIC <- n * sum(obs^2)
    }
    else {
       STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1, n - lag)*obs^2))
    }

    if(sqrd.res) {
       fitdf <- 0;
       names(STATISTIC) <- "X-squared on Squared Residuals for detecting nonlinear processes"
    }
    else if(log.sqrd.res) { 
       fitdf <- 0;
       names(STATISTIC) <- "X-squared on Log-Squared Residuals for detecting nonlinear processes";
    }
    else if(abs.res) { 
       fitdf <- 0;
       names(STATISTIC) <- "X-squared on Absolute valued Residuals for detecting nonlinear processes";
    }
    else {
       names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
    }

    mydf <- lag - fitdf;

    PARAMETER <- c(mydf);
    names(PARAMETER) <- c("df")

    PVAL <- 1 - pchisq(STATISTIC, df=mydf)
    names(PVAL) <- "p-value"

 }
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
       p.value = PVAL, method = METHOD, data.name = DNAME), 
       class = "htest")
}

######
## Weighted.LM.test -- The Weighted Li-Mak type test
##
## Quick intro for the input parameters
##
##      x       -- the residuals or initial unfitted time series
##     h.t      -- the sample conditional variances
##     lag      -- The lag we wish to test at, this is m or M in the literature, default=1
##     type     -- The type of test, here just based on acfs or partial-acfs, default="correlation"
##    fitdf     -- The number of ARCH parameter that have been fit to the series x, default=1
##                 **Note** If fitdf=0, h.t is empty, you want to be performing one of the test in Weighted.Box.test() 
##
## For backward compatability
##
##   weighted   -- For backward compatability, If set to FALSE, you will perform the Li-Mak test, default=TRUE
######

Weighted.LM.test <- function (x, h.t, lag = 1, type = c("correlation", "partial"), fitdf = 1, weighted=TRUE) 
{
 ### Error Checking
 ###
    if (NCOL(x) > 1) 
       stop("x is not a vector or univariate time series");
    if (fitdf >= lag)
       stop("Lag must exceed fitted degrees of freedom");
    if (fitdf < 1)
       stop("Fitted degrees of freedom must be positive");
    if( !(length(x)==length(h.t)) )
       stop("Length of x and h.t must match");

    DNAME <- deparse(substitute(x))
    type <- match.arg(type)

    x <- x^2/h.t;

    if( type == "partial") {
       cor <- acf(x, lag.max = lag, type="partial", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[1:lag];
    }
    else {
       cor <- acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
       obs <- cor$acf[2:(lag + 1)];
    }


    if(type == "correlation" && weighted) {
       METHOD <- "Weighted Li-Mak test on autocorrelations (Gamma Approximation)"
    }
    else if(type == "partial" && weighted) {
       METHOD <- "Weighted Li-Mak test on partial autocorrelations (Gamma Approximation)"
    }
    else if(type == "correlation" && !weighted) {
       METHOD <- "Li-Mak test on autocorrelations (Chi-Squared Approximation)"
    }
    else {
       METHOD <- "Li-Mak test on partial autocorrelations (Chi-Squared Approximation)"
    }

    n <- sum(!is.na(x))
    if(weighted) {
       weights <- (lag - (fitdf+1):lag + (fitdf+1) )/lag;
       obs <- obs[(fitdf+1):lag];
       STATISTIC <- n * sum(weights*obs^2);
       names(STATISTIC) <- "Weighted X-squared on Squared Residuals for fitted ARCH process";
       shape <- (3/4)*(lag + fitdf + 1)^2*(lag - fitdf)/(2*lag^2 + 3*lag + 2*lag*fitdf + 2*fitdf^2 + 3*fitdf + 1);
       scale <- (2/3)*(2*lag^2 + 3*lag + 2*lag*fitdf + 2*fitdf^2 + 3*fitdf + 1)/(lag*(lag + fitdf + 1));
       PARAMETER <- c(shape, scale);
       names(PARAMETER) <- c("Shape", "Scale")
    }
    else {
       weights <- rep(1,(lag-fitdf) );
       obs <- obs[(fitdf+1):lag];
       STATISTIC <- n * sum(weights*obs^2);
       names(STATISTIC) <- "X-squared on Squared Residuals for fitted ARCH process"
       shape <- (lag-fitdf)/2;          # Chi-squared df in Gamma form.
       scale <- 2;
       PARAMETER <- c((lag-fitdf));
       names(PARAMETER) <- c("Degrees of Freedom");
    }
  
    PVAL <- 1 - pgamma(STATISTIC, shape=shape, scale=scale)
    names(PVAL) <- "Approximate p-value"

    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
       p.value = PVAL, method = METHOD, data.name = DNAME), 
       class = "htest")
}
