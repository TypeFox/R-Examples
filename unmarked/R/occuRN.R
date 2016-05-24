
# Fit the Occupancy model of Royle and Nichols

occuRN <- function(formula, data, K = 25, starts, method = "BFGS",
    se = TRUE, ...)
{
    if(!is(data, "unmarkedFrameOccu"))
        stop("Data is not an unmarkedFrameOccu object.")

    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if(is.null(X.offset))
        X.offset <- rep(0, nrow(X))
    if (is.null(V.offset))
        V.offset <- rep(0, nrow(V))

  y <- truncateToBinary(data@y)

  J <- ncol(y)
  M <- nrow(y)

  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)

  nP <- nDP + nOP
  if(!missing(starts) && length(starts) != nP)
	   stop(paste("The number of starting values should be", nP))

  y.ji <- as.vector(y)
  navec <- is.na(y.ji)
  n <- 0:K

  nll <- function(parms, f = "Poisson")
  {

    ## compute individual level detection probabilities
    r.ij <- matrix(plogis(V %*% parms[(nOP + 1) : nP] + V.offset), M, J,
      byrow = TRUE)

    ## compute list of detection probabilities along N
    p.ij.list <- lapply(n, function(k) 1 - (1 - r.ij)^k)

    ## compute P(y_{ij} | N) (cell probabilities) along N
    cp.ij.list <- lapply(p.ij.list, function(pmat) pmat^y * (1-pmat)^(1-y))

    ## replace NA cell probabilities with 1.
    cp.ij.list <- lapply(cp.ij.list, function(cpmat) {
      cpmat[navec] <- 1
      cpmat
    })

    ## multiply across J to get P(y_i | N) along N
    cp.in <- sapply(cp.ij.list, rowProds)

    ## compute P(N = n | lambda_i) along i
    lambda.i <- exp(X %*% parms[1 : nOP] + X.offset)
    lambda.in <- sapply(n, function(x) dpois(x, lambda.i))

    ## integrate over P(y_i | N = n) * P(N = n | lambda_i) wrt n
    like.i <- rowSums(cp.in * lambda.in)

    -sum(log(like.i))
  }

	if(missing(starts)) starts <- rep(0, nP)
  fm <- optim(starts, nll, method = method, hessian = se, ...)
	opt <- fm
	if(se) {
            tryCatch(covMat <- solve(fm$hessian),
                     error=function(x) stop(simpleError("Hessian is singular.  Try providing starting values or using fewer covariates.")))
	} else {
            covMat <- matrix(NA, nP, nP)
	}
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP # + 2 * nP * (nP + 1) / (M - nP - 1)
  names(ests) <- c(occParms, detParms)

  stateEstimates <- unmarkedEstimate(name = "Abundance",
      short.name = "lam",
      estimates = ests[1:nOP],
      covMat = as.matrix(covMat[1:nOP,1:nOP]), invlink = "exp",
      invlinkGrad = "exp")

  detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
      estimates = ests[(nOP + 1) : nP],
      covMat = as.matrix(covMat[(nOP + 1) : nP, (nOP + 1) : nP]),
      invlink = "logistic", invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(state=stateEstimates,
          det=detEstimates))

  umfit <- new("unmarkedFitOccuRN", fitType = "occuRN",
      call = match.call(), formula = formula, data = data,
      sitesRemoved = designMats$removed.sites, estimates = estimateList,
      AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll)

  return(umfit)
}
