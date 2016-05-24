###########################################################################/**
# @RdocFunction fitSnpNmfArray
#
# @title "Allele-specific copy number estimation using non-negative matrix factorization (NMF)"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{Y}{An Lx2xI @array where L is number of probe pairs,
#     2 is the number of alleles (A and B),
#     and I is the number of arrays.}
#  \item{maxIter}{A positive @integer specifying the maximum number of
#     iterations used to calculate the decomposition.}
#  \item{acc}{A positive @double specifying the converence threshold. For
#     more details on convergence, see below.}
# }
#
# \value{
#  Returns a @list of class \code{SnpNmfFit}:
#  \item{Y}{The Lx2xI @array \code{Y}.}
#  \item{W}{The Kx2 @matrix containing allele-specific affinity estimates
#     where K=2L.}
#  \item{H}{A 2xI @matrix containing allele-specific copy number estimates.}
#  \item{hasConverged}{@TRUE if the algorithm converged, otherwise @FALSE.
#     If not applicable, it is @NA.}
#  \item{nbrOfIterations}{The number of iteration ran before stopping.
#     If not applicable, it is @NA.}
# }
#
# \details{
#   The algorithm is considered to have converged when the maximum update
#   of any allele-specific copy number of any array (\code{H}) is greater
#   than \code{acc}.
# }
#
# @examples "../incl/fitSnpNmfArray.Rex"
#
# \seealso{
#   Internally, the array is stacked into a 2LxI matrix and decomposed
#   using @see "fitSnpNmf".
#   See @see "plot.SnpNmfFit".
# }
#
# @keyword internal
#*/###########################################################################
fitSnpNmfArray <- function(Y, ...) {
  # Argument 'Y':
  dim <- dim(Y);
  if (length(dim) != 3L) {
    dimStr <- paste(dim, collapse="x");
    stop("Argument 'Y' is not a three-dimensional array: ", dimStr);
  }

  if (dim[2] != 2L) {
    stop("Second dimension of argument 'Y' is not of length 2: ", dim[2L]);
  }

  # Transform to a "stacked" matrix
  V <- Y;
  dim(V) <- c(dim[1L]*dim[2L], dim[3L]);

  fit <- fitSnpNmf(V, ...);
  fit$V <- V
  fit$Y <- Y
  W2 <- fit$W;
  dim(W2) <- c(dim[1:2], dim(W2)[2L]);
  fit$W2 <- W2;
  fit$args <- list(...);
  class(fit) <- c(class(fit), "SnpNmfFit");

  fit;
} # fitSnpNmfArray()


############################################################################
# HISTORY:
# 2009-03-25 [HB]
# o Added fitSnpNmfArray() accepting a Lx2xI array instead of a 2LxI matrix.
# o Created.
############################################################################
