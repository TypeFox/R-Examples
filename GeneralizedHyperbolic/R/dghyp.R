### Function to calculate the density of the
### generalized hyperbolic distribution
dghyp <- function(x, mu = 0, delta = 1, alpha = 1, beta = 0,
                  lambda = 1, param = c(mu,delta,alpha,beta,lambda)) {

  # Lambda defaults to one if omitted from param vector
  if (length(param) == 4)
    param <- c(param,1)

  ## check parameters
  parResult <- ghypCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  mu <- param[1]
  delta <- param[2]
  alpha <- param[3]
  beta <- param[4]
  lambda <- param[5]
  gamma <- sqrt(alpha^2 - beta^2)

  ## Argument of Bessel K function in numerator
  y <- alpha*sqrt(delta^2 + (x - mu)^2)
  bx <- beta*(x - mu)

  ## Deal with underflow in ratio of Bessel K functions
  ## besselK underflows for x > 740
  ## Use exponentially scaled besselK
  if (delta*gamma > 700) {
    # underflow in constant part
    expTerm <- exp(delta*gamma - y + bx)
    besselRatio <- besselK(x = y, nu = lambda - 1/2, expon.scaled = TRUE)/
                   besselK(x = delta*gamma, nu = lambda, expon.scaled = TRUE)
    expAndBessel <- expTerm*besselRatio
  } else {
    expAndBessel <- ifelse(y > 700 | bx > 700, # underflow in variable part
                           exp(delta*gamma - y + bx)*
                           besselK(x = y, nu = lambda - 1/2,
                                   expon.scaled = TRUE)/
                           besselK(x = delta*gamma, nu = lambda,
                                   expon.scaled = TRUE),
                           exp(bx)*besselK(x = y, nu = lambda - 1/2)/
                           besselK(x = delta*gamma, nu = lambda))
  }

  dens <- (y/alpha)^(lambda - 1/2)*((gamma/delta)^lambda) *
          alpha^(1/2 - lambda)*expAndBessel/sqrt(2*pi)

  dens
} ## End of dghyp()

pghyp <- function (q, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                   param = c(mu,delta,alpha,beta,lambda),
                   lower.tail = TRUE, subdivisions = 100,
                   intTol = .Machine$double.eps^0.25,
                   valueOnly = TRUE, ...)
{
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

   modeDist <- ghypMode(param = param)
   qLess <- which((q <= modeDist)&(is.finite(q)))
   qGreater <- which((q > modeDist)&(is.finite(q)))
   prob <- rep(NA, length(q))
   err <- rep(NA, length(q))

   prob[q == -Inf] <- 0
   prob[q == Inf] <- 0
   err[q %in% c(-Inf, Inf)] <- 0

   dghypInt <- function(q)
   {
       dghyp(q, param = param)
   }

   for (i in qLess)
   {
       intRes <- integrate(dghypInt, -Inf, q[i],
                           subdivisions = subdivisions,
                           rel.tol = intTol, ...)
       prob[i] <- intRes$value
       err[i] <- intRes$abs.error
   }

   for (i in qGreater)
   {
       intRes <- integrate(dghypInt, q[i], Inf,
                           subdivisions = subdivisions,
                           rel.tol = intTol, ...)
       prob[i] <- intRes$value
       err[i] <- intRes$abs.error
   }

   if (lower.tail == TRUE)
   {
       prob[q > modeDist] <- 1 - prob[q > modeDist]
   }
   else
   {
       prob[q <= modeDist] <- 1 - prob[q <= modeDist]
   }

   ifelse(valueOnly, return(prob),
          return(list(value = prob, error = err)))
}

qghyp <- function (p, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                   param = c(mu,delta,alpha,beta,lambda),
                   lower.tail = TRUE, method = c("spline", "integrate"),
                   nInterpol = 501, uniTol = .Machine$double.eps^0.25,
                   subdivisions = 100, intTol = uniTol, ...)
{
    parResult <- ghypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error")
        stop(errMessage)
    if(!lower.tail){
      p <- 1 - p
      lower.tail == TRUE
    }
    method <- match.arg(method)
    mu <- param[1]
    delta <- param[2]
    alpha <- param[3]
    beta <- param[4]
    lambda <- param[5]
    modeDist <- ghypMode(param = param)
    pModeDist<- pghyp(modeDist, param = param, intTol = intTol)
    xRange <- ghypCalcRange(param = param, tol = 10^(-5))

    quant <- rep(NA, length(p))
    invalid <- which((p < 0) | (p > 1))
    pFinite <- which((p > 0) & (p < 1))


    if (method == "integrate"){
        less <- which((p <= pModeDist) & (p > .Machine$double.eps^8))
        quant <- ifelse(p <= .Machine$double.eps^8, -Inf, quant)
        if (length(less) > 0){
            pLow <- min(p[less])
            xLow <- modeDist - sqrt(ghypVar(param = param))
            while (pghyp(xLow, param = param, intTol = intTol) >= pLow)
            {
                xLow <- xLow - sqrt(ghypVar(param = param))
            }
            xRange <- c(xLow, modeDist)
            zeroFn <- function(x, param, p)
            {
                return(pghyp(x, param = param,
                             subdivisions = subdivisions,
                             intTol = intTol) - p)
            }
            for (i in less)
            {
                quant[i] <- uniroot(zeroFn, param = param, p = p[i],
                                    interval = xRange, tol = uniTol)$root
            }
        }

        greater <- which ((p > pModeDist) & (p < (1 - .Machine$double.eps^8)))
        p[greater] <- 1 - p[greater]
        quant <- ifelse(p >=(1 - .Machine$double.eps^8), Inf, quant)
        if (length(greater) > 0){
            pHigh <- min(p[greater])
            xHigh <- modeDist + sqrt(ghypVar(param = param))
            while (pghyp(xHigh, param = param, intTol = intTol,
                         lower.tail = FALSE) >= pHigh)
            {
                xHigh <- xHigh + sqrt(ghypVar(param = param))
            }
            xRange <- c(modeDist, xHigh)
            zeroFn <- function(x, param, p)
            {
                return(pghyp(x, param = param, lower.tail = FALSE,
                             subdivisions = subdivisions,
                             intTol = intTol) - p)
            }
            for (i in greater)
            {
                quant[i] <- uniroot(zeroFn, param = param, p = p[i],
                                    interval = xRange, tol = uniTol)$root
            }
        }
    } else if (method == "spline"){
        inRange <- which((p > pghyp(xRange[1],
                                    param = param, intTol = intTol)) &
                         (p < pghyp(xRange[2], param = param, intTol = intTol)))
        small <- which((p <= pghyp(xRange[1],
                        param = param, intTol = intTol)) & (p >= 0))
        large <- which((p >= pghyp(xRange[2],
                        param = param, intTol = intTol)) & (p <= 1))
        extreme <- c(small, large)
        xVals <- seq(xRange[1], xRange[2], length.out = nInterpol)
        yVals <- pghyp(xVals, param = param, subdivisions = subdivisions,
                       intTol = intTol)
        splineFit <- splinefun(xVals, yVals)
        zeroFn <- function(x, p){
            return(splineFit(x) - p)
        }

        for (i in inRange){
            quant[i] <- uniroot(zeroFn, p = p[i],
                                interval = xRange, tol = uniTol)$root
        }

        if (length(extreme) > 0){
            quant[extreme] <- qghyp(p[extreme], param = param,
                                    lower.tail = lower.tail,
                                    method = "integrate",
                                    nInterpol = nInterpol, uniTol = uniTol,
                                    subdivisions = subdivisions,
                                    intTol = intTol, ...)
        }
    }
    return(quant)
}

### Derivative of the density
ddghyp <- function(x, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                   param = c(mu, delta, alpha, beta, lambda)) {

  # Lambda defaults to one if omitted from param vector
  if (length(param) == 4)
    param <- c(param, 1)

  ## check parameters
  parResult <- ghypCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  mu <- param[1]
  delta <- param[2]
  alpha <- param[3]
  beta <- param[4]
  lambda <- param[5]

  ## Terms for simplification of programming
  t1 <- sqrt(delta^2 + (x - mu)^2)
  t2 <- sqrt(alpha^2 - beta^2)
  t3 <- besselK(x = alpha*t1, nu = lambda - 0.5)
  t4 <- besselK(x = alpha*t1, nu = lambda + 0.5)

  ddghyp <- (t3*(beta*delta^2 + (2*lambda - 1)*(x - mu) + beta*(x - mu)^2) -
             t4*alpha*t1*(x - mu)) *
             exp(beta*(x - mu))*t1^(lambda - (5/2))*t2^lambda/
             (sqrt(2*pi)*alpha^(lambda -  1/2)*delta^lambda *
             besselK(x = delta*t2, nu = lambda))

  ddghyp
} ## End of ddghyp()



### Function to generate random observations from a
### (generalized) hyperbolic distribution using the
### mixing property of the generalized inverse
### Gaussian distribution.
rghyp <- function(n, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                  param = c(mu, delta, alpha, beta, lambda)) {

  # Lambda defaults to one if omitted from param vector
  if (length(param) == 4)
    param <- c(param, 1)

  ## check parameters
  parResult <- ghypCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  mu <- param[1]
  delta <- param[2]
  alpha <- param[3]
  beta <- param[4]
  lambda <- param[5]

  chi <- delta^2
  psi <- alpha^2 - beta^2

  if (lambda == 1) {
    X <- rgig1(n, param = c(chi, psi, lambda))
  } else {
    X <- rgig(n, param = c(chi, psi, lambda))
  }

  sigma <- sqrt(X)
  Z <- rnorm(n)
  Y <- mu + beta*sigma^2 + sigma*Z

  Y
} ## End of rghyp()
