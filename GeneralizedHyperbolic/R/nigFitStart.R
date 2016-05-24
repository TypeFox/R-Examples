### CYD 30/05/10
### DJS 11/09/06
nigFitStart <- function(x, startValues = c("FN","Cauchy","MoM","US"),
                        paramStart = NULL,
                        startMethodMoM = c("Nelder-Mead","BFGS"),
                        ...) {
  startValues <- match.arg(startValues)
  startMethodMoM <- match.arg(startMethodMoM)

  ## right = FALSE deals better with discrete data such as counts
  histData <- hist(x, plot = FALSE, right = FALSE, ...)
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

  if (startValues == "Cauchy") {
    require(MASS)
    svName <- "Cauchy"
    cauchyFit <- fitdistr(x, "Cauchy")
    mu <- cauchyFit$estimate[[1]]
    delta <- cauchyFit$estimate[[2]]
    beta <- 0
    alpha <- 0.1
    paramStart <- c(mu, delta, alpha, beta)
  }


  if (!(startValues %in% c("US", "FN", "Cauchy")))
    startValues <- "MoM"

  if (startValues == "MoM") {
    svName <- "Method of Moments"
    paramStart <- nigFitStartMoM(x, startMethodMoM = startMethodMoM, ...)
  }



  names(paramStart) <- c("mu", "delta", "alpha", "beta")
  list(paramStart = paramStart, breaks = breaks, midpoints = midpoints,
       empDens = empDens, svName = svName)
} ## End of nigFitStart()



nigFitStartMoM <- function(x, startMethodMoM = "Nelder-Mead", ...) {

  fun1 <- function(expParam) {
    diff1 <- nigMean(param = expParam) - mean(x)
    diff1
  }

  fun2 <- function(expParam) {
    diff2 <- nigVar(param = expParam) - var(x)
    diff2
  }

  fun3 <- function(expParam) {
    diff3 <- nigSkew(param = expParam) - skewness(x)
    diff3
  }

  fun4 <- function(expParam) {
    diff4 <- nigKurt(param = expParam) - kurtosis(x)
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
  mu <- mean(x) - delta*hyperbPi*RLambda(zeta, lambda = -1/2)
  startValuesMoM <- c(mu, log(delta), hyperbPi, log(zeta))

  ## Get Method of Moments estimates
  MoMOptim <- optim(startValuesMoM, MoMOptimFun, method = startMethodMoM, ...)
  paramStart <- MoMOptim$par
  paramStart <- hyperbChangePars(1, 2,
                param = c(paramStart[1], exp(paramStart[2]), paramStart[3],
                exp(paramStart[4])))
  return(paramStart)
} ## End of nigFitStartMoM
