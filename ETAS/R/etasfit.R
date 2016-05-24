
etasfit <- function(theta, revents, rpoly, tperiod, integ0, verbose)
{
  tht <- sqrt(theta)
  storage.mode(revents) <- storage.mode(rpoly) <- "double"
  rdata <- list(revents, rpoly, as.double(tperiod), as.double(integ0))
  cfit <- .Call("fit", as.double(tht), rdata, as.integer(verbose), PACKAGE="ETAS")
  list(estimate=cfit[[1]]^2, loglik=cfit[[2]], gradient=cfit[[3]], aic=cfit[[4]])
}

