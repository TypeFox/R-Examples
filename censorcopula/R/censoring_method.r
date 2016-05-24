intervalFitb <- function(copula, method, x, start=NULL, lower=NULL, upper=NULL, 
                       optim.control=list(maxit=1000), 
                       estimate.variance=NA, hideWarnings=TRUE, 
                       bound.eps=.Machine$double.eps^0.5) {
  # use to get fit copula using ml, only log likelihood function is newloglik22
  # Args:
  #   copula: the chosen copula function for simulation
  #   x: round data from the chosen copula function
  #
  # Returns:
  #   estimated parameter of copula function 
  u <- pobs(x)
  # Error handling
  if (any(u < 0) || any(u > 1)) {
    stop("'u' must be in [0,1] -- probably rather use pobs(.)")
  }
  stopifnot(is.numeric(d <- ncol(u)), d >= 2)
  if (copula@dimension != d) {
    stop("The dimension of the data and copula do not match")
  }
  if (is.null(start)) {
    start <- copula:::fitCopStart(copula, u)
  }
  
  if (any(is.na(start))) {
    stop("'start' contains NA values")
  }
  q <- length(copula@parameters)
  if (q != length(start)) {
    stop(gettextf("The lengths of 'start' (= %d) and copula@parameters 
                  (=%d) differ", length(start), q), domain = NA)
  }
  control <- c(optim.control, fnscale = -1)
  control <- control[!vapply(control, is.null, NA)]
  if (!is.null(optim.control[[1]])) {
    control <- c(control, optim.control)
  }
  meth.has.bounds <- method %in% c("Brent", "L-BFGS-B")
  lower <- ifelse (meth.has.bounds, copula@param.lowbnd + bound.eps, -Inf )
  upper <- ifelse (meth.has.bounds, copula@param.upbnd - bound.eps, Inf )
  (if (hideWarnings) { 
    suppressWarnings
  } else {
    identity
  })(fit <- optim(start, Newloglik2, lower = lower, upper = upper,
                  method = method, copula = copula, x = x, control = control))
  return(fit)
  }

Newloglik2 <- function(param, x, copula) {
  # use to get new log likelihood where censoring is taken into consideration
  # Args:
  #   param: parameter of copula
  #   x: round data from the chosen copula function
  #   copula: the chosen copula function for simulation
  copula@parameters <- param
  DupInd <- function(x) {
    # get the duplicated data of x
    # Returns: 
    #   a logical vector indicating which elements of x are duplicates
    return <- duplicated(x) | duplicated(x, fromLast=TRUE)
  }
  maxP <- apply(x, 2, rank, ties.method = "max") / (nrow(x) + 1)
  minP <- apply(x, 2, rank, ties.method = "min") / (nrow(x) + 1)
  Ind <- apply(maxP, 2, DupInd)
  n <- nrow(x)
  cc <- dd <- cd <- dc <- NULL
  # c means the variate at this position has ties
  # e.g. dc: the first X2 has ties
  for (i in 1:n) {
    index <- Ind[i, ]
    maxR <- maxP[i, ]
    minR <- minP[i, ]
    if (!(index[1] | index[2])) {
      dd <- rbind(dd, c(minR, maxR))
    } else if (index[1] & index[2]) {
      cc <- rbind(cc, c(minR, maxR))
    } else if (!index[1]) {
      dc <- rbind(dc, c(minR, maxR))
    } else { 
      cd <- rbind(cd, c(minR, maxR))
    }
  }
  if (!is.null(cc)) {
    CC <- pCopula(cc[,3:4],copula) + pCopula(cc[,1:2],copula) -
      pCopula(cbind(cc[,1], cc[, 4]), copula) - 
      pCopula(cbind(cc[,3], cc[,2]), copula) 
  } else {
    CC <- 1
  }
  if (!is.null(dd)) {
    DD <- dCopula(dd[,1:2], copula) 
  } else {
    DD <- 1
  }
  if (!is.null(cd)) {
    CD <- copula:::dCdu(copula, matrix(cd[, 3:4], ncol=2))[, 2] - 
      copula:::dCdu(copula, matrix(cd[, 1:2], ncol=2))[, 2] 
  } else { 
    CD <- 1
  }
  if (!is.null(dc)) {
    DC <- copula:::dCdu(copula, matrix(dc[, 3:4], ncol=2))[, 1] -
      copula:::dCdu(copula, matrix(dc[, 1:2], ncol=2))[, 1] 
  } else {
    DC <- 1
  }
  return((sum(log(CC)) + sum(log(CD)) + sum(log(DC)) + sum(log(DD))))
}