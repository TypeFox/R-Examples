
##FIXME TODO for other methods (concentrated likelihood, EM)

##FIXME if xreg is not null (in ellipsis, list(...)) and derivatives is "analytical"
#change to numerical (at least for maxlik.fd/td.optim functions)

#stsmFit <- function(x, xreg = NULL, method = c("maxlik.fd.scoring"), args = NULL)
#stsmFit <- function(x, xreg = NULL, method = c("maxlik.fd.scoring"), ...)
stsmFit <- function(x, stsm.method = c("maxlik.fd.scoring", "maxlik.td.scoring", 
  "maxlik.fd.optim", "maxlik.td.optim"), xreg = NULL, ...)
{
  method <- match.arg(stsm.method)

  ##NOTE 
  # argument "xreg" could be avoided because it is supposed to be defined 
  # in the object "x", but this argument is convenient for package "tsoutlies".
  # In this way, the same arguments are passed to "stats::arima" and "stsmFit" 
  # and the code is simplified there without the need of "if" statements.

  if (!is.null(xreg))
  {
    if (!is.null(x@xreg))
      stop("multiple definitions of ", sQuote("xreg"))
    #x@xreg <- xreg
    x <- set.xreg(x, xreg)
  }
  
  #do.call(method, args = c(list(m = x, xreg = xreg), args))
  res <- do.call(stsm.method, args = c(list(m = x), list(...)))
  res$call0 <- res$call
  res$call <- NULL

  # the call to "stsmFit" (not to "maxlik.td.optim" or other "stsm.method")
  # is required in package "tsoutliers::remove.outliers"
  res <- c(call = match.call(), res)
  class(res) <- "stsmFit"

  res
}

# do.call("stsmFit", args = list(x = x, args = list()))
# fcn <- function(x, xreg = NULL, method = c("maxlik.fd.scoring"), ...)
# {
#   #match.call(expand.dots = TRUE)
#   res <- list(...)
#   res
# }
