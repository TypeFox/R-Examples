
decluster <- function(theta, rbwd, revents, rpoly, tperiod)
{
  tht <- sqrt(theta)
  storage.mode(revents) <- storage.mode(rpoly) <- "double"
  cdeclust <- .Call("declust", as.double(tht), as.double(rbwd),
                    revents, rpoly, as.double(tperiod), PACKAGE="ETAS")
  list(revents=cdeclust[[1]], integ0=cdeclust[[2]])
}

