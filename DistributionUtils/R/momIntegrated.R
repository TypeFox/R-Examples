momIntegrated <- function(densFn = "ghyp", param = NULL,
                          order,  about = 0, absolute = FALSE, ...) {

  if (missing(densFn) | !(is.function(densFn) | is.character(densFn)))
    stop("'densFn' must be supplied as a function or name")

  ## Set default integration limits
  low <- -Inf
  high <- Inf

  if (is.character(densFn)) {

    if (is.null(densFn))
      stop("unsupported distribution")
    if (densFn == "ghyp" | densFn == "hyperb" |
        densFn == "gig" | densFn == "vg")
    {
      if (!exists(paste("d", densFn, sep = ""), mode = "function"))
        stop("Relevant package must be loaded")
    }



    if (densFn == "invgamma" | densFn == "inverse gamma"){
      l <- list(...)
      shape <- l$shape
      if(shape <= order)
        stop("Order must be less than shape parameter for inverse gamma")
      low <- 0
      dinvgamma <- function(x, shape, rate = 1, scale = 1/rate) {
        dens <- ifelse(x <= 0, 0,
                       (scale/x)^shape*exp(-scale/x)/(x*gamma(shape)))
        return(dens)
      }

      if (!absolute) {
        ddist <- function(x, order, about, ...) {
          (x - about)^order*dinvgamma(x, ...)
        }
      } else {
        ddist <- function(x, order, about, ...) {
          abs(x - about)^order*dinvgamma(x, ...)
        }
      }
    } else {
      dfun <- match.fun(paste("d", densFn, sep = ""))
      if (densFn == "gamma"){
        l <- list(...)
        shape <- l$shape
        if(order <= -(shape))
          stop("Order must be greater than shape parameter for gamma")
        low <- 0
      }
      if (!absolute) {
        if (is.null(param)){
          ddist <- function(x, order, about, ...) {
            (x - about)^order*dfun(x, ...)
          }
        } else {
          ddist <- function(x, order, about, param) {
            (x - about)^order*dfun(x, param = param)
          }
        }
      } else {
        if (is.null(param)){
          ddist <- function(x, order, about, ...) {
            abs(x - about)^order*dfun(x, ...)
          }
        } else {
          ddist <- function(x, order, about, param) {
            abs(x - about)^order*dfun(x, param = param)
          }
        }
      }
    }
  }
  if (is.null(param)){
    mom <- integrate(ddist, low, high,
                     order = order, about = about,
                     subdivisions = 1000,
                     rel.tol = .Machine$double.eps^0.5, ...)[[1]]
  } else {
    mom <- integrate(ddist, low, high, param = param,
                     order = order, about = about,
                     subdivisions = 1000,
                     rel.tol = .Machine$double.eps^0.5)[[1]]
  }

  ## Return Value:
  return(mom)
}
