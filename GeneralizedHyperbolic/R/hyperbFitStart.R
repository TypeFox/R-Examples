### CYD 01/04/10

hyperbFitStart <- function(x, startValues = c("BN","US","FN","SL","MoM"),
                           paramStart = NULL,
                           startMethodSL = c("Nelder-Mead","BFGS"),
                           startMethodMoM = c("Nelder-Mead","BFGS"),
                           ...) {
  startValues <- match.arg(startValues)
  startMethodSL <- match.arg(startMethodSL)
  startMethodMoM <- match.arg(startMethodMoM)

  ## right = FALSE deals better with discrete data such as counts
  histData <- hist(x, breaks = "Sturges", plot = FALSE, right = FALSE, ...)
  breaks <- histData$breaks
  midpoints <- histData$mids
  empDens <- ifelse(!is.finite(log(histData$density)), NA, histData$density)
  maxIndex <- order(empDens, na.last = FALSE)[length(empDens)]

  if (length(na.omit(empDens[1:maxIndex])) > 1) {
    leftAsymptote <- lm(log(empDens)[1:maxIndex] ~ midpoints[1:maxIndex])$coef
    rightAsymptote <- c(NA, -10*leftAsymptote[2]) # arbitrary large value
  }

  if (length(na.omit(empDens[maxIndex:length(empDens)])) > 1) {
    rightAsymptote <- lm(log(empDens)[maxIndex:length(empDens)] ~
                      midpoints[maxIndex:length(empDens)])$coef

    if (length(na.omit(empDens[1:maxIndex])) < 2)
      leftAsymptote <- c(NA, -10*rightAsymptote[2]) # arbitrary large value
  }

  if ((length(na.omit(empDens[1:maxIndex])) < 2) &
      (length(na.omit(empDens[maxIndex:length(empDens)])) < 2)) {
    if (startValues == "BN" | startValues == "SL")
      stop("not enough breaks to estimate asymptotes to log-density")
  }

  if (startValues == "US") {
    svName <- "User Specified"

    if (is.null(paramStart))
      stop("paramStart must be specified")

    if (!is.null(paramStart)) {
      if (length(paramStart) != 4)
        stop("paramStart must contain 4 values")
      if (paramStart[3] <= 0)
        stop("alpha must be greater than zero")
      if (abs(paramStart[4]) >= paramStart[3])
        stop("absolute value of beta must be less than alpha")
    }
    paramStart <- c(mu = paramStart[1], delta = paramStart[2],
                      alpha = paramStart[3], beta = paramStart[4])
  }

  if (startValues == "FN") {
    svName <- "Fitted Normal"
    nu <- as.numeric(midpoints[maxIndex])
    mu <- mean(x)
    delta <- sd(x)
    hyperbPi <- (nu - mu)/delta
    zeta <- 1 + hyperbPi^2
    paramStart <- hyperbChangePars(1, 2, param = c(mu, delta, hyperbPi, zeta))
  }

  if (startValues == "SL") {
    svName <- "Skew Laplace"

    llsklp <- function(param) {
      -sum(log(dskewlap(x, param = param)))
    }

    lSkewAlpha <- log(1/leftAsymptote[2])
    lSkewBeta <- log(abs(1/rightAsymptote[2]))
    skewMu <- midpoints[maxIndex]
    paramStart <- c(skewMu, lSkewAlpha, lSkewBeta)
    skewlpOptim <- optim(paramStart, llsklp, NULL,
                         method = startMethodSL, hessian = FALSE, ...)
    phi <- 1/exp(skewlpOptim$par[2])
    hyperbGamma <- 1/exp(skewlpOptim$par[3])
    delta <- 0.1 # Take delta to be small
    mu <- skewlpOptim$par[1]
    paramStart <- hyperbChangePars(3, 2, param = c(mu, delta, phi, hyperbGamma))
  }

  if (startValues == "MoM") {
    svName <- "Method of Moments"
    paramStart <- hyperbFitStartMoM(x, startMethodMoM = startMethodMoM, ...)
  }

  if (!(startValues %in% c("US", "FN", "SL", "MoM")))
    startValues <- "BN"

  if (startValues=="BN") {
    svName <- "Barndorff-Nielsen 1977"
    phi <- max(0.001, leftAsymptote[2]) # protect against negative
    hyperbGamma <- max(0.001,-rightAsymptote[2]) # protect against negative

    if (!(is.na(leftAsymptote[1]) | is.na(rightAsymptote[1]))) {
      mu <- -(leftAsymptote[1] - rightAsymptote[1]) /
             (leftAsymptote[2] - rightAsymptote[2])
      intersectionValue <- leftAsymptote[1] + mu*leftAsymptote[2]
      logModalDens <- log(max(empDens, na.rm = TRUE))
      zeta <- intersectionValue - logModalDens

      if (zeta <= 0)
        zeta <- 0.1        # This is set arbitrarily
    } else {
      mu <- median(x)
      intersectionValue <- mu
      logModalDens <- log(max(empDens, na.rm = TRUE))
      zeta <- intersectionValue - logModalDens

      if (zeta <= 0)
        zeta <- 0.1        # This is set arbitrarily
    }

    delta <- zeta/sqrt(phi*hyperbGamma)
    ## Ensure delta is not zero
    if (delta <= 0) {
      delta <- 0.001 # This is set arbitrarily
    }
    paramStart <- hyperbChangePars(3, 2, param = c(mu, delta, phi, hyperbGamma))
  }

  names(paramStart) <- c("mu", "delta", "alpha", "beta")
  list(paramStart = paramStart, breaks = breaks, midpoints = midpoints,
       empDens = empDens, svName = svName)
} ## End of hyperbFitStart()



hyperbFitStartMoM <- function(x, startMethodMoM = "Nelder-Mead", ...) {

  fun1 <- function(expParam) {
    diff1 <- hyperbMean(param = expParam) - mean(x)
    diff1
  }

  fun2 <- function(expParam) {
    diff2 <- hyperbVar(param = expParam) - var(x)
    diff2
  }

  fun3 <- function(expParam) {
    diff3 <- hyperbSkew(param = expParam) - skewness(x)
    diff3
  }

  fun4 <- function(expParam) {
    diff4 <- hyperbKurt(param = expParam) - kurtosis(x)
    diff4
  }

  MoMOptimFun <- function(param) {
    expParam <- hyperbChangePars(1, 2,
                  param = c(param[1], exp(param[2]), param[3], exp(param[4])))
    (fun1(expParam))^2 + (fun2(expParam))^2 +
    (fun3(expParam))^2 + (fun4(expParam))^2
  }

  ## Determine starting values for parameters using
  ## Barndorff-Nielsen et al "The Fascination of Sand" in
  ## A Celebration of Statistics pp.78--79
  xi <- sqrt(kurtosis(x)/3)
  chi <- skewness(x)/3  # Ensure 0 <= |chi| < xi < 1

  if (xi >= 1)
    xi <- 0.999

  adjust <- 0.001

  if (abs(chi) > xi) {
    if (xi < 0 ) {
      chi <- xi + adjust
    } else {
      chi <- xi - adjust
    }
  }

  hyperbPi <- chi/sqrt(xi^2 - chi^2)
  zeta <- 1/xi^2 - 1
  rho <- chi/xi
  delta <- (sqrt(1 + zeta) - 1)*sqrt(1 - rho^2)
  mu <- mean(x) - delta*hyperbPi*RLambda(zeta, lambda = 1)
  startValuesMoM <- c(mu, log(delta), hyperbPi, log(zeta))

  ## Get Method of Moments estimates
  MoMOptim <- optim(startValuesMoM, MoMOptimFun, method = startMethodMoM, ...)
  paramStart <- MoMOptim$par
  paramStart <- hyperbChangePars(1, 2,
                param = c(paramStart[1], exp(paramStart[2]), paramStart[3],
                exp(paramStart[4])))
  return(paramStart)
} ## End of hyperbFitStartMoM
