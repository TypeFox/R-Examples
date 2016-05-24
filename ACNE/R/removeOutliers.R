###########################################################################/**
# @RdocFunction removeOutliers
#
# @title "Removes outliers in matrix containing SNP signals"
#
# \description{
#   @get "title" by identifying outlier elements.  The values of the
#   elements that are outliers are substituted by corresponding values
#   predicted values \code{Yest=W*H} from the current affinity (\code{W})
#   and copy number (\code{H}) estimates.
# }
#
# @synopsis
#
# \arguments{
#  \item{Y}{An IxK @matrix.}
#  \item{W}{A Kx2 @matrix of probe-affinity estimates.}
#  \item{H}{A 2xI @matrix of allele-specific copy-number estimates.}
#  \item{tau}{A scalar specifying the threshold for identifying outliers.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an IxK @matrix where outliers have been "pruned".
#   Outliers are substituted by the corresponding value of \code{Yest}.
# }
#
# @keyword internal
#*/###########################################################################
removeOutliers <- function(Y, W, H, tau=10, ...) {
  # Number of arrays
  I <- ncol(Y);
  # Number of probes
  K <- nrow(Y);

  # Sanity check (may be removed in the future /HB 2009-03-24)
  stopifnot(nrow(W) == K && ncol(W) == 2L);
  stopifnot(nrow(H) == 2L && ncol(H) == I);

  # Output matrix
  Yprime <- Y;

  # Calculating residuals (E) of model Y = W*H + E.
  Yest <- W %*% H;
  E <- Y - Yest;

  # Identify outliers
  rowMad <- rowMads(E);
  outliers <- which(abs(E) > tau*rowMad);

  # Replacing outliers
  if (length(outliers) > 0L) {
    Yprime[outliers] <- Yest[outliers];
  }

  # Sanity check (may be removed in the future /HB 2009-03-24)
  stopifnot(nrow(Yprime) == K && ncol(Yprime) == I);

  Yprime;
} # removeOutliers()


############################################################################
# HISTORY:
# 2009-03-24 [HB]
# o Renamed from RemoveOutliers() to removeOutliers().
# o Added Rdoc comments.
# o Cleaning up code. Minor speed up.
# 2009-02-05 [MO]
# o Fix the value we assign to the outlier
# 2009-02-02 [MO]
# o More efficient code and different name for the indexes (i -> ii)
# 2009-01-30 [MO]
# o Change all the file. Faster and more appropiate for our case
# o Now we use the estimation of the W and the H
# 2009-01-28 [HB]
# o Added an explicit return value.
############################################################################
