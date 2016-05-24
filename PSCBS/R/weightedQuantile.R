############################################################################/**
# @RdocDefault weightedQuantile
#
# @title "Weighted Quantile Value"
#
# @synopsis
#
# \description{
#   Computes a weighted quantile of a numeric vector.
# }
#
# \arguments{
#   \item{x}{a @numeric @vector containing the values whose weighted
#            quantile is to be computed.}
#   \item{w}{a numeric @vector of weights the same length as
#            \code{x} giving the weights to use for each element of \code{x}.
#            Negative weights are treated as zero weights.
#            Default value is equal weight to all values.}
#   \item{probs}{a @numeric @vector of quantiles in [0,1] to be retrieved.}
#   \item{na.rm}{a @logical value indicating whether @NA values in
#            \code{x} should be stripped before the computation proceeds,
#            or not.}
#   \item{method}{If \code{"wtd.quantile"}, then @see "Hmisc::wtd.quantile"
#            of the \pkg{Hmisc} package is used.
#            No other methods are currently supported.}
#   \item{...}{Additional arguments passed to the estimator.}
# }
#
# \value{
#   Returns the weighted quantile.
# }
#
# @author "HB"
#
# \seealso{
#   Internally the following functions may be used:
#   @see "stats::quantile" (if no weights are specified), or
#   @see "Hmisc::wtd.quantile".
#   For a weighted median estimator, @see "matrixStats::weightedMedian"
#   of the \pkg{matrixStats} package.
# }
#
# @keyword univar
# @keyword robust
# @keyword internal
#*/############################################################################
setMethodS3("weightedQuantile", "default", function(x, w, probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE, method=c("wtd.quantile"), ...) {
  # Argument 'x':
  x <- Arguments$getNumerics(x);

  # Argument 'w':
  if (missing(w)) {
    # By default use weights that are one.
    w <- rep(1, times=length(x));
  } else {
    w <- Arguments$getNumerics(w, range=c(0,Inf), length=rep(length(x), times=2L));
  }

  naValue <- NA;
  storage.mode(naValue) <- storage.mode(x);

  # Argument 'na.rm':
  if (is.na(na.rm)) {
    # There are no NAs
  } else if (isTRUE(na.rm)) {
    # Remove values that are NA's
    tmp <- !(is.na(x) | is.na(w));
    x <- .subset(x, tmp);
    w <- .subset(w, tmp);
  } else if (anyNA(x)) {
    return(naValue);
  }

  # Argument 'method':
  method <- match.arg(method);
  if (method == "wtd.quantile") {
    # This will load 'Hmisc', if not already done
    wtd.quantile <- Hmisc::wtd.quantile;
  }



  # Remove values with zero (and negative) weight. This will:
  # (1) take care of the case when all weights are zero,
  # (2) it will most likely speed up the sorting.
  n <- length(w);
  tmp <- (w > 0);
  if (!all(tmp)) {
    x <- .subset(x, tmp);
    w <- .subset(w, tmp);
    n <- sum(tmp);
  }

  # Are there any values left to calculate the weighted median of?
  if (n == 0) {
    return(naValue);
  } else if (n == 1) {
    return(x);
  }

  # Are any weights Inf? Then treat them with equal weight and all others
  # with weight zero. If they have equal weight, regular quantile
  # can be used instead, which is assumed to be faster.
  tmp <- is.infinite(w);
  if (any(tmp)) {
    x <- .subset(x, tmp);
    # Here we know there are no NAs.
    return(quantile(x, probs=probs, na.rm=FALSE, ...));
  }

  # Here we know that there are no missing values in the data
  if (method == "wtd.quantile") {
    wtd.quantile(x, weights=w, probs=probs, normwt=TRUE, na.rm=FALSE, ...);
  } else {
    throw("Cannot estimate weighted quantiles: Argument 'method' is unknown: ", method);
  }
}) # weightedQuantile()


############################################################################
# HISTORY:
# 2013-09-26 [HB]
# o CLEANUP: Now weightedQuantile(..., method=="wtd.quantile") no longer
#   attaches 'Hmisc', but only loads its namespace.
# 2012-08-30
# o Updated Rdoc cross reference for matrixStats to point to matrixStats.
# 2011-04-08
# o Created.
############################################################################
