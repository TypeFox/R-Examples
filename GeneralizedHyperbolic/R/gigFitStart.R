### Function to find start values for fitting the GIG distribution
### David Cusack, 4/6/2010
### David Scott 2/11/2010

gigFitStart <- function(x, startValues = c("LM","GammaIG","MoM","Symb","US"),
                        paramStart = NULL,
                        startMethodMoM = c("Nelder-Mead","BFGS"), ...) {

  startValues <- match.arg(startValues)
  startMethodMoM <- match.arg(startMethodMoM)

  histData <- hist(x, plot = FALSE, right = FALSE, ...)
  breaks <- histData$breaks
  midpoints <- histData$mids
  empDens <- ifelse(!is.finite(log(histData$density)), NA, histData$density)

  if (startValues == "LM"){
    svName <- "Linear Models"
    paramStart <- gigFitStartLM(x, ...)
  }

  if (startValues == "GammaIG") {
    svName <- "Gamma and Inverse Gamma"
    lg <- try(fitdistr(x, "gamma")$loglik, silent = TRUE)
    lig <- try(fitdistr(1/x, "gamma")$loglik, silent = TRUE)

    c <- class(lg)
    d <- class(lig)

    if (c == "try-error") {
      if (d == "try-error") {
        stop("GammaIG method does not work for this data set")
      }   else g <- fitdistr(1/x, "gamma")$estimate
      paramStart <- c(2*g[1], 0.1, -1/g[2])
    } else if (d == "try-error") {
      if (c == "try-error") {
        stop("GammaIG method does not work for this data set")
      }   else g <- fitdistr(x, "gamma")$estimate
      paramStart <- c(0.1, 2*g[1], g[2])
    } else if (lg > lig) {
      g <- fitdistr(x, "gamma")$estimate
      paramStart <- c(0.1, 2*g[1], g[2])
    } else if (lg < lig) {
      g <- fitdistr(1/x, "gamma")$estimate
      paramStart <- c(2*g[1], 0.1, -1/g[2])
    }
  }

  if (startValues == "MoM") {
    svName <- "Method of Moments"
    paramStart <- gigFitStartMoM(x, startMethodMoM = startMethodMoM, ...)
  }

  if (startValues == "Sym") {
    svName <- "Symbolic"
    stop("Symbolic method not yet implemented")
  }

  if (startValues == "US")  {
    svName <- "User Specified"
    if (is.null(paramStart))
      stop("paramStart must be specified")
    if (!is.null(paramStart)) {
      if (length(paramStart) != 3)
        stop("paramStart must contain 3 values")
      if (paramStart[1] < 0)
        stop("chi must be greater than 0")
      if (paramStart[2] < 0)
        stop("psi must be greater than 0")
    }
  }

  names(paramStart) <- c("chi","psi","lambda")
  list(paramStart = paramStart, breaks = breaks, midpoints = midpoints,
       empDens = empDens, svName = svName)
}

### Method of moments
gigFitStartMoM <- function(x, paramStart = NULL,
                           startMethodMoM = "Nelder-Mead", ...) {

  ## Define functions to fit moments
  fun1 <- function(expParam) {
    diff1 <- mean(x) - gigMean(param = expParam)
    diff1
  }
  fun2 <- function(expParam) {
    diff2 <- var(x) - gigVar(param = expParam)
    diff2
  }
  fun3 <- function(expParam) {
    diff3 <- skewness(x) - gigSkew(param = expParam)
    diff3
  }

  MoMOptimFunc <- function(param) {
    expParam <- c(exp(param[1]), exp(param[2]), param[3])
    return((fun1(expParam))^2 + (fun2(expParam))^2 + (fun3(expParam))^2)
    }

  ## Fitting log(chi), log(psi) and lambda
  if (is.null(paramStart)){
    startValuesMoM <- c(1,1,1)
  }
  paramStart <- c(log(paramStart[1]),log(paramStart[2]),paramStart[3])

  MoMOptim <- optim(startValuesMoM, MoMOptimFunc,
                    method = startMethodMoM, ...)
  paramStart <- c(exp(MoMOptim$par[1]), exp(MoMOptim$par[2]),
                     MoMOptim$par[3])
  ## Small chi and psi can cause problems later
  if (paramStartMoM[1] < 0.1)
    paramStartMoM <- c(0.1, paramStartMoM[2], paramStartMoM[3])
  if (paramStartMoM[2] < 0.1)
        paramStartMoM <- c(paramStartMoM[1], 0.1, paramStartMoM[3])
  paramStartMoM

  return(paramStart)
}

### Using linear models fits to tails
gigFitStartLM <- function(x, ...)
{
  ## Purpose: Find starting values for fitting the gig distribution
  ## ----------------------------------------------------------------------
  ## Arguments: x, data
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 28 Oct 2010, 14:07

  ## Call to histogram
  histInfo <- hist(x, plot = FALSE, ...)
  mids <- as.matrix(histInfo$mids)

  ## Find density
  dens <- density(x)

  ## Default values
  lambda1 <- NA
  lambda2 <- NA
  chi <- NA
  psi <- NA

  ## Find the closest density estimate to the midpoint
  func <- function(x) which.min(abs(dens$x - x))
  estx <- apply(mids, 1, func)
  esty <- dens$y[estx]
  nMax <- which.max(esty)

  ## Only keep cells not below the mode & counts not too small
  nCells <- length(mids)
  aboveMode <- rep(TRUE, nCells)
  aboveMode[1:(nMax - 1)] <- FALSE
  keep <- rep(FALSE, nCells)
  keep <- ifelse((log(esty) > -9)&aboveMode, TRUE, keep)

  lx <- log(estx)[keep]
  invx <- (1/estx)[keep]
  ly <- log(esty)[keep]

  ## Linear model
  if (length(keep[keep==TRUE]) >=2){
    estPars <- coef(lm(ly ~ estx[keep]  + lx ))[-1]
    psi <- -2*estPars[1]
    lambda1 <- estPars[2] + 1
  }

  ## Now invert x
  x <- 1/x

  histInfo <- hist(x, plot = FALSE, ...)
  mids <- as.matrix(histInfo$mids)

  ### Find density
  dens <- density(x)

  ## Find the closest density estimate to the midpoint
  func <- function(x) which.min(abs(dens$x - x))
  estx <- apply(mids, 1, func)
  esty <- dens$y[estx]
  nMax <- which.max(esty)

  ## Only keep cells not below the mode & counts not too small
  nCells <- length(mids)
  aboveMode <- rep(TRUE, nCells)
  aboveMode[1:(nMax - 1)] <- FALSE
  keep <- rep(FALSE, nCells)
  keep <- ifelse((log(esty) > -9)&aboveMode, TRUE, keep)

  lx <- log(estx)[keep]
  invx <- (1/estx)[keep]
  ly <- log(esty)[keep]

  ## Linear model
  if (length(keep[keep==TRUE])>= 2){
    estPars <- coef(lm(ly ~ estx[keep]  + lx ))[-1]
    chi <- -2*estPars[1]
    lambda2 <- -estPars[2] - 1
  }
  lambda <- (lambda1 + lambda2)/2
  ### Deal with special cases
  if (is.na(lambda1)) lambda <- lambda2
  if (is.na(lambda2)) lambda <- lambda1
  if (is.na(lambda1)&is.na(lambda2)) lambda <- 1

  ## Ensure parameters are not outside range
  if (is.na(chi)) chi <- 1
  if (is.na(psi)) psi <- 1
  if (chi < (.Machine$double.eps)^0.5) chi <- 0.01
  if (psi < (.Machine$double.eps)^0.5) psi <- 0.01
  paramStart <- c(chi,psi,lambda)
  return(paramStart)
}
