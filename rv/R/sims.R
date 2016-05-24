

sims <- function(x, ...) {
  UseMethod("sims")
}

sims.rvsummary <- function(x, dimensions=FALSE,  ...) {
  ## NOEXPORT
  S <- matrix(unlist(x, use.names=FALSE), nrow=length(x[[1]]))
  q <- attr(x, "quantiles")
  if (!is.null(q)) {
    q.names <- paste(100*q, "%", sep="")
    dimnames(S) <- list(q.names, names(x))
  } else if (!is.null(a <- attr(x, "levels"))) {
    dimnames(S) <- list(a, names(x))
  }
  if (dimensions) {
    if (!is.null(dim. <- dim(x))) {
      S <- t(S)
      dim(S) <- c(dim., dim(S)[2])
      ndim <- length(dim(S))
      newdim <- c(ndim, 1:(ndim-1))
      S <- aperm(S, newdim)
    }
  }
  attr(S, "quantiles") <- q
  return(S)
}

sims.default <- function(x, ...) {
  ## NOEXPORT
  ##
  ## We should not set the dimension of constants.
  ## re: problems in max.rvsims: cbind(sims(1),sims(x)) does not work
  ## since sims(1) returns a 1x1 matrix but sims(x) may be L x n matrix.
  ## TODO: integrate better!
  ## Need to unclass first, to be consistent.
  ## DEBUG: ? why not just sims(as.rv(x, ...)) ? 
  if (is.atomic(x)) {
    as.vector(unclass(x))  # drop attributes
  } else {
    stop("No method to extract the simulations of an object of class ", paste(class(x), collapse="/"))
  }
}

# ========================================================================
# sims.rv  -  return the matrix of simulations for an rv
# ========================================================================
# gives the simulations of a vector as a matrix,
# used e.g. for computing quantiles of a vector.
# length(x) columns, rvnsims rows.
# length of sims is adjusted to the longest sims length
# if dimensions=TRUE, will return an L x dim(x) matrix.
# (this is useful for simapply applied to random matrices, e.g. see determinant.rv)
#
#
 ###, as.list=FALSE)



sims.rv <- function(x, dimensions=FALSE, n.sims=getnsims(), ...) {
  ## NOEXPORT
  if (length(x)<1) {
    return(matrix(nrow=n.sims,ncol=0))
  }
  xl <- rvnsims(x)
  if (missing(n.sims)) {
    n.sims <- max(xl)
  } else if (!is.numeric(n.sims) && n.sims>=2) {
    stop("n.sims (if specified) must be at least 2")
  }
  class(x) <- NULL
  for (i in which(xl!=n.sims)) {
    x[[i]] <- rep(x[[i]], length.out=n.sims)
  }
  m <- matrix(unlist(x), nrow=n.sims)
  if (dimensions && !is.null(dim.x <- dim(x))) {
    n.s <- (length(m)%/%prod(dim.x))
    cm <- m
    dimnames(cm) <- NULL
    m <- array(m, c(n.s, dim.x)) # multi-way array, 1st dimension is the dimension "rows"
    if (!is.null(dn <- dimnames(x))) {
      dimnames(m)[1+seq_along(dn)] <- dn ## put dimnames in correct positions
    } 
  } else {
    dimnames(m) <- list(NULL, names(x))
  }
  return(m)
}

