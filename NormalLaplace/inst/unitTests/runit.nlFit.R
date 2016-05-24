### Testing nlFit
test.nlFit <- function() {
  param <- c(0, 1, 2, 3)
  names(param) <- c("mu", "sigma", "alpha", "beta")

  ## RUnit uses kind = "M-M", normal.kind = "K-R" for RNG. See ?RNGkind
  set.seed(2242, kind = "Marsaglia-Multicarry")
  dataVector <- rnl(1000, param = param)

  ## Grabbing starting parameter object
  testParamStart <- nlFitStart(dataVector)

  ## Running fitting over different methods of optimisation
  testnlFitDefault <- nlFit(dataVector, hessian = TRUE)
  testnlFitBFGS <- nlFit(dataVector, method = "BFGS", hessian = TRUE)
  testnlFitLBFGSB <- nlFit(dataVector, method = "L-BFGS-B", hessian = TRUE)
  testnlFitnlm <- nlFit(dataVector, method = "nlm", hessian = TRUE)
  testnlFitnlminb <- nlFit(dataVector, method = "nlminb", hessian = TRUE)

  ## These checks will either pass or fail regardless of method choice
  checkEquals(dataVector, testnlFitDefault$obs)
  checkEquals("dataVector", testnlFitDefault$obsName)
  checkEquals(is.null(testnlFitDefault$hessian), FALSE)
  checkEquals(testParamStart$paramStart, testnlFitDefault$paramStart)
  checkEquals(testParamStart$svName, testnlFitDefault$svName)
  checkEquals(testParamStart$breaks, testnlFitDefault$breaks)
  checkEquals(testParamStart$empDens, testnlFitDefault$empDens)
  checkEquals(testParamStart$midpoints, testnlFitDefault$midpoints)
  checkEquals(testParamStart$startValues, testnlFitDefault$startValues)

  ## Checking default/N-M
  checkEquals("Nelder-Mead", testnlFitDefault$method)
  checkTrue(is.numeric(testnlFitDefault$conv))
  checkTrue(is.numeric(testnlFitDefault$iter) &
            testnlFitDefault$iter > 0)  # More than one iteration occurred
  checkTrue(is.numeric(testnlFitDefault$maxLik))
  checkTrue(all(! is.na(testnlFitDefault$param)) &
            is.numeric(testnlFitDefault$param))  # Want no NAs and all numeric
  checkEquals(class(testnlFitDefault), c("nlFit", "distFit"))

  ## Checking BFGS
  checkEquals("BFGS", testnlFitBFGS$method)
  checkTrue(is.numeric(testnlFitBFGS$conv))
  checkTrue(is.numeric(testnlFitBFGS$iter) &
            testnlFitBFGS$iter > 0)  # More than one iteration occurred
  checkTrue(is.numeric(testnlFitBFGS$maxLik))
  checkTrue(all(! is.na(testnlFitBFGS$param)) &
            is.numeric(testnlFitBFGS$param))  # Want no NAs and all numeric
  checkEquals(class(testnlFitBFGS), c("nlFit", "distFit"))

  ## Checking L-BFGS-B
  checkEquals("L-BFGS-B", testnlFitLBFGSB$method)
  checkTrue(is.numeric(testnlFitLBFGSB$conv))
  checkTrue(is.numeric(testnlFitLBFGSB$iter) &
            testnlFitLBFGSB$iter > 0)  # More than one iteration occurred
  checkTrue(is.numeric(testnlFitLBFGSB$maxLik))
  checkTrue(all(! is.na(testnlFitLBFGSB$param)) &
            is.numeric(testnlFitLBFGSB$param))  # Want no NAs and all numeric
  checkEquals(class(testnlFitLBFGSB), c("nlFit", "distFit"))

  ## Checking nlm
  checkEquals("nlm", testnlFitnlm$method)
  checkTrue(is.numeric(testnlFitnlm$conv))
  checkTrue(is.numeric(testnlFitnlm$iter) &
            testnlFitnlm$iter > 0)  # More than one iteration occurred
  checkTrue(is.numeric(testnlFitnlm$maxLik))
  checkTrue(all(! is.na(testnlFitnlm$param)) &
            is.numeric(testnlFitnlm$param))  # Want no NAs and all numeric
  checkEquals(class(testnlFitnlm), c("nlFit", "distFit"))

  ## Checking nlminb
  checkEquals("nlminb", testnlFitnlminb$method)
  checkTrue(is.numeric(testnlFitnlminb$conv))
  checkTrue(is.numeric(testnlFitnlminb$iter) &
            testnlFitnlminb$iter > 0)  # More than one iteration occurred
  checkTrue(is.numeric(testnlFitnlminb$maxLik))
  checkTrue(all(! is.na(testnlFitnlminb$param)) &
            is.numeric(testnlFitnlminb$param))  # Want no NAs and all numeric
  checkEquals(class(testnlFitnlminb), c("nlFit", "distFit"))
}


## Testing graphical output
graphicstest.nlFit <- function() {
  param <- c(0, 1, 2, 3)
  names(param) <- c("mu", "sigma", "alpha", "beta")

  ## RUnit uses kind = "M-M", normal.kind = "K-R" for RNG. See ?RNGkind
  set.seed(2242, kind = "Marsaglia-Multicarry")
  dataVector <- rnl(1000, param = param)
  testnlFit <- nlFit(dataVector)

  pdf("Histogram of dataVector.pdf")
  plot.nlFit(testnlFit, which = 1)
  dev.off()

  pdf("Log-Histogram of dataVector.pdf")
  plot.nlFit(testnlFit, which = 2)
  dev.off()
}
