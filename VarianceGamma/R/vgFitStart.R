vgFitStart <- function (x, breaks = NULL, startValues = "SL",
                        paramStart = NULL, startMethodSL = "Nelder-Mead",
                        startMethodMoM = "Nelder-Mead", ...) {
  histData <- hist(x, plot = FALSE, right = FALSE)
    if (is.null(breaks)) {
      breaks <- histData$breaks
    }
    midpoints <- histData$mids
    empDens <- ifelse(!is.finite(log(histData$density)), NA, histData$density)
    maxIndex <- order(empDens, na.last = FALSE)[length(empDens)]
    if (length(na.omit(empDens[1:maxIndex])) > 1) {
      leftAsymptote <- lm(log(empDens)[1:maxIndex] ~
                          midpoints[1:maxIndex])$coef
      rightAsymptote <- c(NA, -10 * leftAsymptote[2])
    }
    if (length(na.omit(empDens[maxIndex:length(empDens)])) > 1) {
      rightAsymptote <- lm(log(empDens)[maxIndex:length(empDens)] ~
                           midpoints[maxIndex:length(empDens)])$coef
      if (length(na.omit(empDens[1:maxIndex])) < 2) {
        leftAsymptote <- c(NA, -10 * rightAsymptote[2])
      }
    }
    if ((length(na.omit(empDens[1:maxIndex])) < 2) &
      (length(na.omit(empDens[maxIndex:length(empDens)])) < 2)) {
      if (startValues == "SL")
        stop("not enough breaks to estimate asymptotes to log-density")
    }

    if (startValues == "US") {
        svName <- "User Specified"
        if (is.null(paramStart))
            stop("ThetaStart must be specified")
        if (!is.null(paramStart)) {
            if (length(paramStart) != 4)
                stop("ThetaStart must contain 4 values")
            if (paramStart[2] <= 0)
                stop("sigma in ThetaStart must be greater than zero")
            if (paramStart[4] <= 0)
                stop("nu in ThetaStart must be greater than zero")
        }
        paramStart <- c(vgC = paramStart[1], sigma = log(paramStart[2]),
                        theta = paramStart[3], nu = log(paramStart[4]))
    }

    if (startValues == "SL") {
        svName <- "Skew Laplace"
        llsklp <- function(param) {
            -sum(log(dskewlap(x, param = param, logPars = TRUE)))
        }
        lSkewAlpha <- log(1/leftAsymptote[2])
        lSkewBeta <- log(abs(1/rightAsymptote[2]))
        skewMu <- midpoints[maxIndex]
        paramStart <- c(lSkewAlpha, lSkewBeta, skewMu)
        skewlpOptim <- optim(paramStart, llsklp, NULL, method = startMethodSL,
                             hessian = FALSE, ...)
        phi <- 1/exp(skewlpOptim$par[1])
        hyperbGamma <- 1/exp(skewlpOptim$par[2])
        delta <- 0.1
        mu <- skewlpOptim$par[3]
        hyperbAlpha <- hyperbChangePars(3, 2,
                                        c(mu, delta, phi, hyperbGamma))[1]
        hyperbBeta <- hyperbChangePars(3, 2,
                                       c(mu, delta, phi, hyperbGamma))[2]
        nu <- 1
        vgC <- mu
        squareSigma <- 2*(1/nu)/(hyperbAlpha^2 - hyperbBeta^2)
        if (squareSigma < 0){
          sigma <- 0.001
        } else {
            sigma <- sqrt(2*(1/nu)/(hyperbAlpha^2 - hyperbBeta^2))
        }
        theta <- hyperbBeta*squareSigma
        paramStart <- c(vgC, log(sigma), theta, log(nu))
    }

    if (startValues == "MoM") {
        svName <- "Method of Moments"
        paramStart <- vgFitStartMoM(x, startMethodMoM = startMethodMoM, ...)
    }

    names(paramStart) <- c("vgC", "lsigma", "theta", "lnu")
    list(paramStart = paramStart, breaks = breaks, midpoints = midpoints,
         empDens = empDens, svName = svName)
}

vgFitStartMoM <- function (x, startMethodMoM = "Nelder-Mead", ...) {
    fun1 <- function(expParam) {
      diff1 <- vgMean(param = expParam) - mean(x)
      return(diff1)
    }
    fun2 <- function(expParam) {
      diff2 <- vgVar(param = expParam) - var(x)
      return(diff2)
    }
    fun3 <- function(expParam) {
      diff3 <- vgSkew(param = expParam) - skewness(x)
      return(diff3)
    }
    fun4 <- function(expParam) {
      diff4 <- vgKurt(param = expParam) - kurtosis(x)
      return(diff4)
    }
    MoMOptimFun <- function(param) {
      expParam <- c(param[1], exp(param[2]), param[3], exp(param[4]))
      return((fun1(expParam))^2 + (fun2(expParam))^2 +
             (fun3(expParam))^2 + (fun4(expParam))^2)
    }

    sigma <- sqrt(var(x))
    nu <- (kurtosis(x) - 3)/3
    adjust <- 0.001
    if (nu < 0 | (abs(nu - 0) < adjust)) {
      nuTrue <- nu    #nu true is the negative value of nu given by the
                      #previous formula if no adjust has been made
      nu <- adjust
      theta <- (skewness(x)*sigma)/(3*nuTrue)
    } else {
      theta <- (skewness(x)*sigma)/(3*nu)
    }
    vgC <- mean(x) - theta

    startValuesMoM <- c(vgC,log(sigma),theta,log(nu))
    MoMOptim <- optim(startValuesMoM, MoMOptimFun,
                      method = startMethodMoM, ...)
    paramStart <- MoMOptim$par
    return(paramStart)
}
