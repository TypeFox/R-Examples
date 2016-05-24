
##NOTE
#"oeff" could be returned as a "ts" but it would require 
# passing that information to the function, for the time being, if needed
# define it as a "ts" outside this function

##NOTE
# argument "pars" is required only for "IOeffect";
# "pars" can be NULL, thus it is not appropriate to use it to dispatch 
# between method "outliers.effects.ArimaPars" or "outliers.effects.stsmSS" 
# as done for example with function "outliers.regressors";
# nevertheless, the argument "tsmethod" can be avoided by using 
# the name of the class of "pars" ("ArimaPars" or "stsmSS")

outliers.effects <- function(mo, n, weights = FALSE, delta = 0.7, 
  pars = NULL, n.start = 50, freq = 12)
{
  # effects of outliers on the original series
  # variables in equation (19) in Chen-Liu (1993)
  # see equations (2) and (4) in the document attached to the package

  IOeffect <- function(n, ind, pars, w = 1, n.start = 50)
  {
    if (!identical(w, 1))
      stopifnot(length(ind) == length(w))

    msg <- paste("argument", sQuote("pars"), "should be of class", 
      sQuote("ArimaPars"), "or", sQuote("stsmSS")) 
      #"\n       see", sQuote("?coefs2poly"), "in package", sQuote("tsoutliers"),
      #"or", sQuote("?char2numeric"), "in package", sQuote("stsm"))

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

      IOeffect(n = n, ind = ind, pars = pars, w = w, n.start = n.start)
    }

    if (inherits(pars, "ArimaPars"))
    {
      mas <- ARMAtoMA(ar = pars$arcoefs, ma = pars$macoefs, 
        lag.max = n - 1)

      m <- matrix(0, nrow = n, ncol = length(ind))
      m[n * seq.int(0, ncol(m) - 1) + ind] <- w
##FIXME see na.omit in the presence of outliers in "x"
      m <- apply(m, 2, function(x, f, nm1) na.omit(filter(c(rep(0, nm1), x), 
        filter = f, method = "conv", sides = 1)), f = c(1, -mas), nm1 = n - 1)
    } else 
    if (inherits(pars, "stsmSS"))
    {
##FIXME is this always unit impulse? if so redo without running KF

##FIXME if for w=1, they are just lagged versions, if so, redo running KF only once
#check it also in tsmethod == "arima"
      m <- matrix(0, nrow = n, ncol = length(ind))
      #m[n * seq.int(0, ncol(m) - 1) + ind] <- w
      #if (w[i] != 1)
      #  f <- f * w[i]
      I <- rep(0, n + n.start)
      for (i in seq_along(ind))
      {
        #I <- ts(rep(0, n + n.start), start = start(resid), frequency = frequency(resid))
        I[] <- 0  
        I[n.start] <- 1
        tmp <- KFKSDS::KF(I, pars)
        f <- tmp$v[-seq.int(n.start-1)]
##FIXME see how to avoid remove last element this way (see doc)
        m[,i] <- f[-length(f)]
      }
    } else
      stop(msg)
    m
  }

  AOeffect <- function(n, ind, w = 1)
  {
    if (!identical(w, 1))
      stopifnot(length(ind) == length(w))
    m <- matrix(0, nrow = n, ncol = length(ind))
    m[n * seq.int(0, ncol(m) - 1) + ind] <- w
    m
  }

  LSeffect <- function(n, ind, w = 1)
  {
    if (!identical(w, 1))
      stopifnot(length(ind) == length(w))
    if (!identical(w, 1)) {
      stopifnot(length(ind) == length(w))
    } else
      w <- rep(1, length(ind))
    m <- matrix(0, nrow = n, ncol = length(ind))
    #m[n * seq.int(0, ncol(m) - 1) + ind] <- w
    #m <- apply(m, 2, function(x) diffinv(x)[-1])
    #m <- apply(m, 2, filter, filter = 1, method = "rec")
    for (i in seq_along(ind))
      m[seq.int(ind[i], n),i] <- w[i]
    m
  }

  TCeffect <- function(n, ind, w = 1, delta = 0.7)
  {
    if (!identical(w, 1))
      stopifnot(length(ind) == length(w))
    m <- matrix(0, nrow = n, ncol = length(ind))
    m[n * seq.int(0, ncol(m) - 1) + ind] <- w
    m <- apply(m, 2, filter, filter = delta, method = "rec")
    m
  }

  SLSeffect <- function(n, ind, freq, w = 1)
  {
    if (!identical(w, 1)) {
      stopifnot(length(ind) == length(w))
    } else
      w <- rep(1, length(ind))
    if (freq == 1 || missing(freq))
      stop(sQuote("freq"), "is either not defined or is not higher than 1")
    m <- matrix(0, nrow = n, ncol = length(ind))
    #m[n * seq.int(0, ncol(m) - 1) + ind] <- w
    #m <- apply(m, 2, filter, filter = c(rep(0, freq - 1), 1), method = "rec")
    for (i in seq_along(ind))
      m[seq.int(ind[i], n, by = freq),i] <- w[i]
    m
  }

  if (!weights)
    mo[,"coefhat"] <- 1

  # "oeff" is filled following the order in which
  # the outliers are defined by rows in 'mo'

  oeff <- matrix(nrow = n, ncol = nrow(mo))

  if (length(indio <- which(mo[,"type"] == "IO")) > 0)
  {
    oeff[,indio] <- IOeffect(n = n, ind = mo[indio,"ind"], pars = pars, 
        w = mo[indio,"coefhat"], n.start = n.start)
  }

  if (length(indao <- which(mo[,"type"] == "AO")) > 0)
  {
    oeff[,indao] <- AOeffect(n, mo[indao,"ind"], mo[indao,"coefhat"])
  }

  if (length(indls <- which(mo[,"type"] == "LS")) > 0)
  {
    oeff[,indls] <- LSeffect(n, mo[indls,"ind"], mo[indls,"coefhat"])
  }

  if (length(indtc <- which(mo[,"type"] == "TC")) > 0)
  {
    oeff[,indtc] <- TCeffect(n, mo[indtc,"ind"], mo[indtc,"coefhat"], delta)
  }

  if (length(indsls <- which(mo[,"type"] == "SLS")) > 0)
  {
    oeff[,indsls] <- SLSeffect(n, mo[indsls,"ind"], freq, mo[indsls,"coefhat"])
  }

  colnames(oeff) <- paste(as.character(mo[,"type"]), mo[,"ind"], sep = "")

  oeff
}
