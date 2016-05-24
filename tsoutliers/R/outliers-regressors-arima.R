
outliers.regressors <- function(pars, mo, n, weights = TRUE,
  delta = 0.7, freq = 12, n.start = 50)
  UseMethod("outliers.regressors")

outliers.regressors.default <- function(pars, mo, n, weights = TRUE,
  delta = 0.7, freq = 12, n.start = 50)
{
  msg <- paste("argument", sQuote("pars"), "should be of class", 
    sQuote("ArimaPars"), "or", sQuote("stsmSS"), 
    "\n       see", sQuote("?coefs2poly"), "in package", sQuote("tsoutliers"),
    "or", sQuote("?char2numeric"), "in package", sQuote("stsm"))

  #if (!inherits(pars, "ArimaPars") && !inherits(pars, "stsmSS"))
  #  stop(msg)

  # if "pars" is a list, deduce the type of time series model 
  # looking at the names of the elements in "pars"

  if (inherits(pars, "list"))
  {
    tmp <- names(pars)
    if (all(c("arcoefs", "macoefs") %in% tmp)) {
      pars <- structure(pars, class = "ArimaPars")
    } else 
    if (all(c("Z", "T", "H") %in% tmp)) {
      pars <- structure(pars, class = "stsmSS")
    } else
      stop(msg)

    outliers.regressors(pars = pars, mo = mo, n = n, weights = weights,
      delta = delta, freq = freq, n.start = n.start)

  } else
    stop(msg)
}

outliers.regressors.ArimaPars <- function(pars, mo, n, weights = TRUE,
  delta = 0.7, freq = 12, n.start = 50)
{
  # "n.start" is not used here with "ArimaPars"
  
  # regressor variables in the auxiliar regression 
  # where outliers are detected
  # equation (20) in Chen-Liu (1993)
  # equation (3) in documentation attached to the package

  # the regressions are not actually run since the t-statistics
  # can be obtained more conveniently as indicated in 
  # equation (14) in Chen-Liu (1993), as it is done in functions 
  # "IOtstat.arima", "AOtstat.arima", "LStstat.arima", 
  # "TCtstat.arima", "SLStstat.arima";
  # these variables are used in function "locate.outliers.iloop" to 
  # adjust the residuals at each iteration

  IOxreg <- function(n, ind, w)
  {
    mI <- matrix(0, nrow = n, ncol = length(ind))
    mI[n * seq.int(0, ncol(mI) - 1) + ind] <- w
    mI
  }

  AOxreg <- function(n, ind, w)
  {
    mI <- matrix(0, nrow = n, ncol = length(ind))
    mI[n * seq.int(0, ncol(mI) - 1) + ind] <- 1

    for (i in seq_along(ind))
    {
      mI[,i] <- w[i] * na.omit(filter(c(rep(0, n-1), mI[,i]), filter = fma, 
        method = "conv", sides = 1))
    }
    mI
  }

  LSxreg <- function(n, ind, w)
  {
    mI <- matrix(0, nrow = n, ncol = length(ind))
    mI[n * seq.int(0, ncol(mI) - 1) + ind] <- 1

    j <- 1
    for (i in ind)
    {
      mI[,j] <- w[j] * diffinv(na.omit(filter(c(rep(0, n-1), mI[,j]), 
        filter = fma, method = "conv", sides = 1)))[-1]
      j <- j + 1
    }
    mI
  }

  TCxreg <- function(n, ind, w, ar, ma, delta = 0.7)
  {   
    mI <- matrix(0, nrow = n, ncol = length(ind))
    mI[n * seq.int(0, ncol(mI) - 1) + ind] <- 1

    madmL <- coef(polynomial(c(1, ma)) * polynomial(c(1, -delta)))[-1]
    ftc <- c(1, ARMAtoMA(-madmL, -ar, n-1))
    
    for (i in seq_along(ind))
    {
      mI[,i] <- w[i] * na.omit(filter(c(rep(0, n-1), mI[,i]), 
        filter = ftc, method = "conv", sides = 1))
    }
    mI
  }

  SLSxreg <- function(n, ind, w, freq)
  {
    mI <- matrix(0, nrow = n, ncol = length(ind))
    mI[n * seq.int(0, ncol(mI) - 1) + ind] <- 1

    j <- 1
    for (i in ind)
    {
      mI[,j] <- w[j] * diffinv(na.omit(filter(c(rep(0, n-1), mI[,j]), 
        filter = fma, method = "conv", sides = 1)), lag = freq)[-seq_len(freq)]
      j <- j + 1
    }
    mI
  }

  if (!weights)
    mo[,"coefhat"] <- 1

  arcoefs <- pars$arcoefs
  macoefs <- pars$macoefs

  if (any(mo[,"type"] %in% c("AO", "LS", "SLS")))
    fma <- c(1, ARMAtoMA(-macoefs, -arcoefs, n-1))

  oxreg <- NULL

  if (length(indio <- which(mo[,"type"] == "IO")) > 0)
  {
    oxreg <- cbind(oxreg, 
      IOxreg(n, mo[indio,"ind"], mo[indio,"coefhat"]))
  }

  if (length(indao <- which(mo[,"type"] == "AO")) > 0)
  {
    oxreg <- cbind(oxreg, 
      AOxreg(n, mo[indao,"ind"], mo[indao,"coefhat"]))
  }

  if (length(indls <- which(mo[,"type"] == "LS")) > 0)
  {
    oxreg <- cbind(oxreg, 
      LSxreg(n, mo[indls,"ind"], mo[indls,"coefhat"]))
  }

  if (length(indtc <- which(mo[,"type"] == "TC")) > 0)
  {
    oxreg <- cbind(oxreg, 
      TCxreg(n, mo[indtc,"ind"], mo[indtc,"coefhat"], arcoefs, macoefs, delta))
  }

  if (length(indsls <- which(mo[,"type"] == "SLS")) > 0)
  {
    oxreg <- cbind(oxreg, 
      SLSxreg(n, mo[indsls,"ind"], mo[indsls,"coefhat"], freq))
  }

  oxreg
}
