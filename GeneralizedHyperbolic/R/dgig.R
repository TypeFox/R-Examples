### Function to calculate the density of the
### generalized inverse Gaussian distribution
dgig <- function(x, chi = 1, psi = 1, lambda = 1,
                 param = c(chi, psi, lambda), KOmega = NULL) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]

  omega <- sqrt(chi * psi)

  if (is.null(KOmega))
    KOmega <- besselK(omega, nu = lambda)

  gigDensity <- ifelse(x > 0, (psi / chi)^(lambda / 2) /
                       (2 * KOmega) * x^(lambda - 1) *
                       exp(-(1/2) * (chi * x^(-1) + psi * x)), 0)
  gigDensity
} ## end of dgig()

### Cumulative distribution function of the generalized inverse Gaussian
### Uses incomplete Bessel function of Slevinsky and Safouhi
###
### DJS 9/8/2010
pgig <- function(q, chi = 1, psi = 1, lambda = 1,
                    param = c(chi,psi,lambda), lower.tail = TRUE,
                    ibfTol = .Machine$double.eps^(0.85),
                    nmax = 200) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]
  qCalculate <- which((q > 0) & (is.finite(q)))
  prob <- rep(NA, length(q))
  ## probabilities are of upper tail: change later if needs be
  prob[q <= 0] <- 1
  prob[q == Inf] <- 0
  omega <- sqrt(chi * psi)
  KOmega <- besselK(omega, nu = lambda)
  x <- psi*q/2
  y <- chi/(2*q)
  const <- (psi/chi)^(lambda/2)/(2*KOmega)

  for (i in qCalculate){
    prob[i] <- q[i]^lambda*incompleteBesselK(x[i], y[i], -lambda,
                                             tol = ibfTol, nmax = nmax)
  }
  prob[qCalculate] <- const*prob[qCalculate]

  if (lower.tail) prob <- 1 - prob

  return(prob)
} ## End of pgig()

### qgig using pgig based on incomplete Bessel function
### David Scott 09/08/2010
qgig <- function(p, chi = 1, psi = 1, lambda = 1,
                 param = c(chi, psi, lambda),
                 lower.tail = TRUE, method = c("spline", "integrate"),
                 nInterpol = 501, uniTol = 10^(-7),
                 ibfTol = .Machine$double.eps^(0.85), nmax =200, ...){

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  if(!lower.tail){
      p <- 1 - p
      lower.tail == TRUE
  }

  method <- match.arg(method)
  param <- as.numeric(param)
  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]
  modeDist <- gigMode(param = param)
  pModeDist<- pgig(modeDist, param = param, ibfTol = ibfTol)
  xRange <- gigCalcRange(param = param, tol = 10^(-7))

  quant <- rep(NA, length(p))
  invalid <- which((p < 0) | (p > 1))
  pFinite <- which((p > 0) & (p < 1))


  if (method == "integrate")
  {
    less <- which((p <= pModeDist) & (p > .Machine$double.eps^5))
    quant <- ifelse(p <= .Machine$double.eps^5, 0, quant)
    if (length(less) > 0){
      ## pLow <- min(p[less])
      ## xLow <- modeDist - sqrt(gigVar(param = param))
      ## while (pgig(xLow, param = param, ibfTol = 10^(-5)) >= pLow){
      ##   xLow <- xLow - sqrt(gigVar(param = param))
      ## }
      xLow <- 0
      xRange <- c(xLow, modeDist)
      zeroFn <- function(x, param, p)
      {
        return(pgig(x, param = param, ibfTol = ibfTol) - p)
      }
      for (i in less){
        quant[i] <- uniroot(zeroFn, param = param, p = p[i],
                            interval = xRange, tol = uniTol)$root
      }
    }

    greater <- which ((p > pModeDist) & (p < (1 - .Machine$double.eps^5)))
    p[greater] <- 1 - p[greater]
    quant <- ifelse(p >= (1 - .Machine$double.eps), Inf, quant)
    if (length(greater) > 0){
      pHigh <- min(p[greater])
      xHigh <- modeDist + sqrt(gigVar(param = param))
      while (pgig(xHigh, param = param, lower.tail = FALSE) >= pHigh){
        xHigh <- xHigh + sqrt(gigVar(param = param))
      }
      xRange <- c(modeDist,xHigh)
      zeroFn <- function(x, param, p)
      {
        return(pgig(x, param = param, lower.tail = FALSE,
                     ibfTol = ibfTol) - p)
      }
      for (i in greater){
        quant[i] <- uniroot(zeroFn, param = param, p = p[i],
                            interval = xRange, tol = uniTol)$root
      }
    }
  } else if (method == "spline") {
    inRange <- which((p > pgig(xRange[1], param = param)) &
                     (p < pgig(xRange[2], param = param)))
    small <- which((p <= pgig(xRange[1], param = param)) & (p >= 0))
    large <- which((p >= pgig(xRange[2], param = param)) & (p <= 1))
    extreme <- c(small, large)
    xVals <- seq(xRange[1], xRange[2], length.out = nInterpol)
    yVals <- pgig(xVals, param = param,
                  ibfTol = max(ibfTol, .Machine$double.eps^0.25) )
    splineFit <- splinefun(xVals, yVals)
    zeroFn <- function(x, p){
      return(splineFit(x) - p)
    }

    for (i in inRange){
      quant[i] <- uniroot(zeroFn, p = p[i],
                          interval = xRange, tol = uniTol)$root
    }

    if (length(extreme) > 0){
      quant[extreme] <- qgig(p[extreme], param = param,
                             lower.tail = lower.tail, method = "integrate",
                             nInterpol = nInterpol, uniTol = uniTol,
                             ibfTol = ibfTol, ...)
    }
  }
  return(quant)
}


# Modified version of rgig to generate random observations
# from a generalized inverse Gaussian distribution in the
# special case where lambda = 1.
rgig1 <- function(n, chi = 1, psi = 1, param = c(chi, psi)) {

  if (length(param) == 2)
    param <- c(param, 1)

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  chi <- param[1]
  psi <- param[2]
  lambda <- 1

  alpha <- sqrt(psi / chi)
  beta <- sqrt(psi * chi)

  m <- abs(beta) / beta

  g <- function(y) {
    0.5 * beta * y^3 - y^2 * (0.5 * beta * m + lambda + 1) +
      y * (-0.5 * beta) + 0.5 * beta * m
  }

  upper <- m

  while (g(upper) <= 0)
    upper <- 2 * upper

  yM <- uniroot(g, interval = c(0, m))$root

  yP <- uniroot(g, interval = c(m, upper))$root

  a <- (yP - m) * exp(-0.25 * beta * (yP + 1 / yP - m - 1 / m))
  b <- (yM - m) * exp(-0.25 * beta * (yM + 1 / yM - m - 1 / m))
  c <- -0.25 * beta * (m + 1 / m)

  output <- numeric(n)

  for (i in 1:n) {
    needValue <- TRUE

    while (needValue) {
      R1 <- runif(1)
      R2 <- runif(1)
      Y <- m + a * R2 / R1 + b * (1 - R2) / R1
      if (Y > 0) {
        if (-log(R1) >= 0.25 * beta * (Y + 1 / Y) + c) {
          needValue <- FALSE
        }
      }
    }

    output[i] <- Y
  }

  output / alpha
} ## End of rgig1

# Function to generate random observations from a
# generalized inverse Gaussian distribution. The
# algorithm is based on that given by Dagpunar (1989)
rgig <- function(n, chi = 1, psi = 1, lambda = 1,
                 param = c(chi, psi, lambda)) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]

  if (lambda == 1)
    stop(return(rgig1(n, param = c(chi, psi))))

  alpha <- sqrt(psi / chi)
  beta <- sqrt(psi * chi)

  m <- (lambda - 1 + sqrt((lambda - 1)^2 + beta^2)) / beta

  g <- function(y) {
    0.5 * beta * y^3 - y^2 * (0.5 * beta * m + lambda + 1) +
      y * ((lambda - 1) * m - 0.5 * beta) + 0.5 * beta * m
  }

  upper <- m

  while (g(upper) <= 0)
    upper <- 2 * upper

  ## yM <- uniroot(g, interval = c(0, m))$root
  ## Correct problem when psi and chi are very small
  ## Code from Fabian Scheipl
  ## Fabian.Scheipl@stat.uni-muenchen.de
  yM <- uniroot(g, interval = c(0, m),
                tol = min(.Machine$double.eps^0.25,
                          (.Machine$double.eps + g(0) / 10)))$root

  yP <- uniroot(g, interval = c(m, upper))$root

  a <- (yP - m) * (yP / m)^(0.5 * (lambda - 1)) *
         exp(-0.25 * beta * (yP + 1 / yP - m - 1 / m))
  b <- (yM - m) * (yM / m)^(0.5 * (lambda - 1)) *
         exp(-0.25 * beta * (yM + 1 / yM - m - 1 / m))
  c <- -0.25 * beta * (m + 1 / m) + 0.5 * (lambda - 1) * log(m)

  output <- numeric(n)

  for (i in 1:n) {
    needValue <- TRUE

    while (needValue) {
      R1 <- runif(1)
      R2 <- runif(1)
      Y <- m + a * R2 / R1 + b * (1 - R2) / R1
      if (Y > 0) {
        if (-log(R1) >= -0.5 * (lambda - 1) * log(Y) +
            0.25 * beta * (Y + 1 / Y) + c) {
          needValue <- FALSE
        }
      }
    }
    output[i] <- Y
  }
  output / alpha
} ## End of rgig()

### Derivative of dgig
ddgig <- function(x, chi = 1, psi = 1, lambda = 1,
                  param = c(chi, psi, lambda), KOmega = NULL) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]

  omega <- sqrt(chi * psi)

  if (is.null(KOmega))
    KOmega <- besselK(x = omega, nu = lambda)

  ddgig <- ifelse(x > 0,
                  dgig(x, param = param, KOmega)*
                  (chi/x^2 + 2*(lambda - 1)/x - psi) / 2,
                  0)
  ddgig
} ## End of ddgig()

