##' Randomly sample response patterns given a list of items
##'
##' Returns a random sample of response patterns given a list of item
##' models and parameters. If \code{grp} is given then theta, items, params,
##' mean, and cov can be omitted.
##'
##' @name rpf.sample
##' @param theta either a vector (for 1 dimension) or a matrix (for >1
##' dimension) of person abilities or the number of response patterns
##' to generate randomly
##' @param items a list of item models
##' @param params a list or matrix of item parameters. If omitted, random item
##' parameters are generated for each item model.
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param prefix Column names are taken from param or items.
##' If no column names are available, some will be generated using
##' the given prefix.
##' @param mean mean vector of latent distribution (optional)
##' @param cov covariance matrix of latent distribution (optional)
##' @param mcar proportion of generated data to set to NA (missing completely at random)
##' @param grp a list with spec, param, mean, and cov
##' @return Returns a data frame of response patterns
##' @export
##' @examples
##' # 1 dimensional items
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' i2 <- rpf.nrm(outcomes=3)
##' i2.p <- rpf.rparam(i2)
##' rpf.sample(5, list(i1,i2), list(i1.p, i2.p))
##' @seealso \code{\link{sample}}
rpf.sample <- function(theta, items, params, ..., prefix="i",
                       mean=NULL, cov=NULL, mcar=0.0, grp=NULL)
{
    garbageArguments <- list(...)
    if (length(garbageArguments) > 0) {
        stop("rpf.sample does not accept values for the '...' argument")
    }

    if (!missing(grp)) {
	    if (missing(items)) items <- grp$spec
	    if (missing(params) && !is.null(grp$param)) params <- grp$param
	    if (missing(mean) && !is.null(grp$mean)) mean <- grp$mean
	    if (missing(cov) && !is.null(grp$cov)) cov <- grp$cov
	    if (missing(theta) && !is.null(grp$data)) theta <- nrow(grp$data)
    }

  numItems <- length(items)
  maxDim <- max(vapply(items, function(i) i@factors, 0))

    if (maxDim > 1) {
	    design <- matrix(rep(1:maxDim, numItems), nrow=maxDim)
	    design[sapply(items, function(i) 1:maxDim > i@factors)] <- NA
    } else {
	    design <- matrix(rep(1, numItems), nrow=1)
    }

  maxAbilities <- max(design, na.rm=TRUE)

  numPeople <- NA
  if (is.numeric(theta) && length(theta) == 1) {
    if (theta <= 1) stop("Request at least 2 samples")
    numPeople <- theta
    if (missing(mean)) mean <- rep(0, maxAbilities)
    if (length(mean) != maxAbilities) stop(paste("Mean vector must have length",maxAbilities,"not",length(mean)))
    if (missing(cov)) cov <- diag(maxAbilities)
    if (any(dim(cov) != maxAbilities)) stop(paste("Cov matrix must be square matrices of size",maxAbilities))
    theta <- array(t(mvtnorm::rmvnorm(numPeople, mean=mean, sigma=cov)),
                   dim=c(maxAbilities, numPeople))
  } else if (maxDim == 1 && is.vector(theta)) {
    numPeople <- length(theta)
    theta <- array(theta, dim=c(maxAbilities, numPeople))
  } else {
      if (dim(theta)[1] > maxAbilities) {
          stop(paste("Only", maxAbilities, "abilities but theta provides", dim(theta)[1],
                     "-- maybe transpose theta?"))
      }
    numPeople <- dim(theta)[2]
  }

  if (missing(params)) {
    params <- lapply(items, rpf.rparam)
  }

  outcomes <- vapply(items, function(i) i@outcomes, 0)
  
  name <- colnames(params)
  if (is.null(name)) name <- names(items)

    maxParam <- max(sapply(items, rpf.numParam))
    if (is.matrix(params) && nrow(params) > maxParam) {
	    warning(paste("Item parameter matrix has", nrow(params) - maxParam, "extra rows"))
    }

  ret <- list()
  for (ix in 1:numItems) {
    i <- items[[ix]]
    param <- c()
    if (is.list(params)) {
      param <- params[[ix]]
    } else {
      param <- params[,ix]  # item parameters are in columns
    }

    if (length(param) < rpf.numParam(i)) {
      stop(paste("Item",class(i),"needs",rpf.numParam(i),
                 "parameters but only",length(param),"given"))
    }

    cols <- design[,ix]
    cols <- cols[!is.na(cols)]
    i.theta <- theta[cols,,drop=FALSE]
    P <- rpf.prob(i, param[1:rpf.numParam(i)], i.theta)
#    if (any(is.na(P))) stop(paste("Item", i@spec, "with param", param," produced NAs"))
    ret1 <- apply(P, 2, sample, x=1:i@outcomes, size=1, replace=F)
    labels <- 1:i@outcomes
    if (!missing(grp) && !is.null(grp$data)) labels <- levels(grp$data[1, name[ix]])
    ret1 <- factor(ret1, levels=1:i@outcomes, ordered=TRUE, labels=labels)
    attr(ret1, 'mxFactor') <- TRUE  # for OpenMx
    ret[[ix]] <- ret1
  }
  ret <- as.data.frame(ret, optional=TRUE)
  if (is.null(name)) name <- paste(prefix,1:numItems,sep="")
  colnames(ret) <- name
  if (mcar > 0) {
      size <- prod(dim(ret))
      mask <- rep(FALSE, size)
      mask[sample.int(size, size * mcar)] <- TRUE
      shaped.mask <- array(mask, dim=dim(ret))
      ret[shaped.mask] <- NA
  }
  return(ret)
}
