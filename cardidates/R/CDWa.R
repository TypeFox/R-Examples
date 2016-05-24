`CDWa` <-
function(p, xmin=0, xmax=365, quantile = 0.05, symmetric = FALSE) {
  xminuser <- xmin
  xmaxuser <- xmax
  # extract parameters if p is a whole cardiFit object
  if (inherits(p, "cardiFit")) p <- p$p
    
  npar <- length(p)
  if (npar == 6) {
    pp       <- c(1, p[2:3], 0, p[5:6]) ## set baseline for integration to zero
    aweibull <- aweibull6
    fweibull <- fweibull6
    xmin <- p[2]  # p[2] = left turnpoint, attention: overwrites user input
    xmax <- p[5]  # p[5] = right turnpoint, attention: overwrites user input
  } else if (npar ==4) {
    pp       <- c(0, p[2:4])
    aweibull <- aweibull4
    fweibull <- fweibull4
  } else stop("invalid number of parameters")

  fzero <- function(x, q, p, x0=0) {
    q - aweibull(lower=x0, upper=x, p)
  }

  tMid <- optimize(f = fweibull, p=p, lower = xmin, upper=xmax, maximum=TRUE)$maximum
  ymax <- fweibull(tMid, p)

  if (symmetric) {
    q        <- c(quantile/2, 1 - quantile/2)
    integral <- aweibull(lower=0, upper = Inf, p=pp)
    tBegin   <- uniroot(fzero, interval=c(0, xmaxuser), q = q[1] * integral, p = pp)$root
    tEnd     <- uniroot(fzero, interval=c(0, xmaxuser), q = q[2] * integral, p = pp)$root
  } else {
    q      <- c(quantile, 1 - quantile)
    left   <- aweibull(lower = 0,    upper = tMid, p=pp)
    right  <- aweibull(lower = tMid, upper = xmaxuser,  p=pp)
    tBegin <- uniroot(fzero, interval=c(0, tMid),   q = q[1] * left,  p = pp)$root
    tEnd   <- uniroot(fzero, interval=c(tMid, xmaxuser), q = q[2] * right, p = pp, x0 = tMid)$root
  }
  list(x = c(tMid=tMid, tBegin=tBegin, tEnd=tEnd),
       y  = fweibull(c(tMid, tBegin, tEnd), p),
       p  = p
  )
}

