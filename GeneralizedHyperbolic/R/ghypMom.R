### calculate moments (mu moments, central moments, raw moments or others)
### for the generalized hyperbolic distribution
### if only momType is defined, then calculate moments according to momType.
### if only about is defined, then calculate moment according to about.
### if users define both momType and about(contradictive or not), then
### the about value would always overwrites momType. Thus, calculate moments
### about "about".

ghypMom <- function(order, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                    param = c(mu, delta, alpha, beta, lambda),
                    momType = c("raw", "central", "mu"), about = 0) {

  ## check order is whole number
  if (!is.wholenumber(order))
    stop("Order must be a whole number")

  if ((order < 0))
    stop("Order must be positive")

  ## check momType
  momType <- match.arg(momType)

  if (! momType %in% c("raw", "central", "mu"))
    stop("Unrecognised moment type")

  ## unpack parameters
  param <- as.numeric(param)

  if (length(param) == 4)
    param <- c(param, 1)

  ## check parameters
  parResult <- ghypCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  mu <- param[1]
  delta <- param[2]
  alpha <- param[3]
  beta <- param[4]
  lambda <- param[5]

  gamma <- sqrt(alpha^2 - beta^2)
  zeta <- delta * gamma

  if (order == 0) {
    mom <- 1
  } else {
    ## calculate mu moments
    muMom <- rep (NA, order)

    for (i in 1:order) {
      a <- momRecursion(order = i)
      coeff <- a$a
      betaPow <- a$M
      deltaPow <- 2 * a$L
      zetaPow <- a$L
      lengthZetaPow <- length(zetaPow)

      ## calculate terms and sum
      muM <- coeff * (delta^deltaPow) * (beta^betaPow) *
        sapply(zetaPow, besselRatio, x = zeta, nu = lambda) / (zeta^zetaPow)
      muMom[i] <- sum(muM)
    }
  }

  if (about != 0) {
    mom <- momChangeAbout(order = order, oldMom = muMom,
                          oldAbout = mu, newAbout = about)
  } else {
    if (momType == "mu")
      mom <- muMom[order]

    if (momType == "raw") {
      about <- 0
      mom <- momChangeAbout(order = order, oldMom = muMom,
                            oldAbout = mu, newAbout = about)
    }

    if (momType == "central") {
      about <- ghypMean(param = param)
      mom <- momChangeAbout(order = order, oldMom = muMom,
                            oldAbout = mu, newAbout = about)
    }
  }

  return(mom)
}
