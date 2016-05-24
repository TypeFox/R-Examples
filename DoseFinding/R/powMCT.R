## all design related functions for power calculations

mvtnorm.control <- function(maxpts = 30000, abseps = 0.001,
                            releps = 0, interval = NULL){
  res <- list(maxpts = maxpts, abseps = abseps,
              releps = releps, interval = interval)
  class(res) <- "GenzBretz"    
  res
}

powCalc <- function(alternative, critV, df, corMat, deltaMat, control){
  nC <- nrow(corMat) # number of contrasts
  if(alternative[1] == "two.sided"){
    lower <- rep(-critV, nC)
  } else {
    lower <- rep(-Inf, nC)
  }
  upper <- rep(critV, nC)
  if (!missing(control)) {
    if(!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  ctrl$interval <- NULL      # not used with pmvt
  nScen <- ncol(deltaMat)
  res <- numeric(nScen)
  for(i in 1:nScen){
    pmvtCall <- c(list(lower, upper, df = df, corr = corMat, delta = deltaMat[,i],
                       algorithm = ctrl))
    res[i] <- as.vector(1 - do.call("pmvt", pmvtCall))
  }
  names(res) <- colnames(deltaMat)
  res
}

powMCT <- function(contMat, alpha = 0.025, altModels,
                   n, sigma, S, placAdj = FALSE,
                   alternative = c("one.sided", "two.sided"),
                   df, critV = TRUE,
                   control = mvtnorm.control()){
  alternative <- match.arg(alternative)
  if(inherits(contMat, "optContr")){
    if(attr(contMat, "placAdj") != placAdj){
      message("using \"placAdj\" specification from contMat object")
      placAdj <- attr(contMat, "placAdj")
    }
    contMat <- contMat$contMat
  }
  if(!is.matrix(contMat))
    stop("contMat needs to be a matrix")
  nD <- nrow(contMat) # nr of doses
  nC <- ncol(contMat) # nr of contrasts
  ## extract covariance matrix
  if(missing(S)){
    if(missing(n) | missing(sigma))
      stop("Either S or n and sigma need to be specified")
    if(length(n) == 1)
      n <- rep(n, nD)
    if(length(n) != nD)
      stop("n needs to be of length nrow(contMat)")
    S <- sigma^2*diag(1/n)
    df <- sum(n) - nD
  } else {
    if(!missing(n)|!missing(sigma))
      stop("Need to specify exactly one of \"S\" or \"n\" and \"sigma\"")
    if(nrow(S) != ncol(S))
      stop("S needs to be a square matrix")
    if(nrow(S) != nD)
      stop("S needs to have as many rows&cols as there are doses")
    if(missing(df))
      stop("need to specify degrees of freedom in \"df\", when specifying \"S\"")
  }
  ## extract means under the alternative
  if(missing(altModels))
    stop("altModels argument needs to be specified")
  muMat <- getResp(altModels)
  if(placAdj){
    muMat <- sweep(muMat, 2, muMat[1,], "-") # remove placebo column
    muMat <- muMat[-1, , drop=FALSE]
  }
  if(nrow(muMat) != nD)
    stop("Incompatible contMat and muMat")
  ## extract df
  if(missing(S)){
    if(missing(df))
      stop("degrees of freedom need to be specified in df")
    df <- sum(n) - nD
  }
  ## calculate non-centrality parameter
  deltaMat <- t(contMat) %*% muMat
  covMat <- t(contMat) %*% S %*% contMat
  den <- sqrt(diag(covMat))
  deltaMat <- deltaMat/den
  if(alternative == "two.sided"){
    deltaMat <- abs(deltaMat)
  }
  corMat <- cov2cor(covMat)
  
  if(!is.finite(df))
    df <- 0
  ## calculate critical value
  if(is.logical(critV) & critV == TRUE){
    critV <- critVal(corMat, alpha, df, alternative, control)
  } # else assume critV already contains critical value
  res <- powCalc(alternative, critV, df, corMat, deltaMat, control)
  ## class(res) <-  "powMCT"
  ## attr(res, "type") <- ifelse(missing(n), "S", "n&sigma")
  ## attr(res, "contMat") <- contMat
  ## attr(res, "muMat") <- muMat
  res
}

## print.powMCT <- function(x, ...){
##   attributes(x)[2:5] <- NULL
##   print(x)
## }
