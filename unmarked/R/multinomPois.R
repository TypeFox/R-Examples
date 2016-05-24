
## Take an M x J matrix of detection probabilities and return a matrix
## of M x J observation probs
# Compute the cell probabilities for the observation classes
# in removal sampling.
#
# Both p and the returned matrix are M x J for M sites and J sampling occasions.

removalPiFun <- function(p){
  M <- nrow(p)
  J <- ncol(p)
  pi <- matrix(NA, M, J)
  pi[,1] <- p[,1]
  for(i in seq(from = 2, length = J - 1)) {
	  pi[, i] <- pi[,i-1] / p[,i-1] * (1-p[,i-1]) * p[,i]
  }
  return(pi)
}

# p is an M x 2 matrix of detection probabilities (site x observer).
# returns an M x 3 matrix of row=(1 not 2, 2 not 1, 1 and 2).
# Compute the cell probabilities for the observation classes
# in double observer sampling.

doublePiFun <- function(p){
  M <- nrow(p)
  pi <- matrix(NA, M, 3)
  pi[,1] <- p[,1] * (1 - p[,2])
  pi[,2] <- p[,2] * (1 - p[,1])
  pi[,3] <- p[,1] * p[,2]
  return(pi)
}


# Fit the multinomial-Poisson abundance mixture model.

multinomPois <- function(formula, data, starts, method = "BFGS",
    se = TRUE, ...)
{
    if(!is(data, "unmarkedFrameMPois"))
		    stop("Data is not a data frame or unmarkedFrame.")
    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    if (is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
        }
    J <- ncol(y)
    R <- obsNum(data)
    M <- nrow(y)
    piFun <- data@piFun

    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)
    nP <- nDP + nAP
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    yvec <- as.numeric(y)
    navec <- is.na(yvec)

    nll <- function(parms) {
        lambda <- exp(X %*% parms[1 : nAP] + X.offset)
        p <- plogis(V %*% parms[(nAP + 1) : nP] + V.offset)
        p.matrix <- matrix(p, M, R, byrow = TRUE)
        pi <- do.call(piFun, list(p = p.matrix))
        logLikeSite <- dpois(y, matrix(lambda, M, J) * pi, log = TRUE)
        logLikeSite[navec] <- 0
        -sum(logLikeSite)
        }
    if(missing(starts))
        starts <- rep(0, nP)
    fm <- optim(starts, nll, method = method, hessian = se, ...)
    opt <- fm
    if(se) {
        tryCatch(covMat <- solve(fm$hessian),
                 error=function(x) stop(simpleError("Hessian is singular.  Try providing starting values or using fewer covariates.")))
    } else {
        covMat <- matrix(NA, nP, nP)
      }
    ests <- fm$par
    fmAIC <- 2 * fm$value + 2 * nP
    names(ests) <- c(lamParms, detParms)

    stateName <- "Abundance"

    stateEstimates <- unmarkedEstimate(name = stateName,
                                       short.name = "lambda",
                                       estimates = ests[1:nAP],
                                       covMat = as.matrix(
                                           covMat[1:nAP,1:nAP]),
                                       invlink = "exp",
                                       invlinkGrad = "exp")

    detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
                                     estimates = ests[(nAP + 1) : nP],
                                     covMat = as.matrix(covMat[(nAP + 1) :
                                                 nP, (nAP + 1) : nP]),
                                     invlink = "logistic",
                                     invlinkGrad = "logistic.grad")

    estimateList <- unmarkedEstimateList(list(state=stateEstimates,
        det=detEstimates))

    umfit <- new("unmarkedFitMPois", fitType = "multinomPois",
        call = match.call(), formula = formula, data = data,
        estimates = estimateList, sitesRemoved = designMats$removed.sites,
        AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll)

    return(umfit)
}
