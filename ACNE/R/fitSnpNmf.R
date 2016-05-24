###########################################################################/**
# @RdocFunction fitSnpNmf
#
# @title "Non-negative matrix factorization (NMF) of a matrix containing SNP probe signals"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{V}{An KxI @matrix where I is the number of arrays and K is the
#     number of probe where K should be even (K=2L).}
#  \item{acc}{A positive @double specifying the converence threshold. For
#     more details on convergence, see below.}
#  \item{maxIter}{A positive @integer specifying the maximum number of
#     iterations used to calculate the decomposition.}
#  \item{maxIterRlm}{A positive @integer specifying the maximum number of
#     iterations used in rlm.}
#  \item{refs}{An index @vector (@integer or @logical) specifying the
#     reference samples.  If @NULL, all samples are used as a reference.}
# }
#
# \value{
#  Returns a @list:
#  \item{W}{The Kx2 @matrix containing allele-specific affinity estimates.}
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
# \seealso{
#   @see "WHInit", @see "robustWInit", @see "robustHInit", and
#   @see "removeOutliers".
# }
#
# @keyword internal
#*/###########################################################################
fitSnpNmf <- function(V, acc=0.02, maxIter=10, maxIterRlm=20, refs=NULL) {
  I <- ncol(V);
  K <- nrow(V);

  # Argument 'refs':
  if (is.null(refs)) {
    refs <- seq(length=I);
  } else if (!is.vector(refs)) {
    throw("Argument 'refs' is not a vector: ", class(refs)[1L]);
  } else if (is.logical(refs)) {
    if (length(refs) != I) {
      throw("The number of elements in argument 'refs' does not match the number of column in argument 'V': ", length(refs), " != ", I);
    }
    refs <- which(refs);
  } else if (is.numeric(refs)) {
    if (!all(1L <= refs & refs <= I)) {
      throw("Some elements in argument 'refs' is out of range [1,", I, "].");
    }
    refs <- as.integer(refs);
  } else {
    throw("Argument 'refs' must be either a logical or a numeric vector: ", mode(refs));
  }

  # A small positive value
  eps <- 1e-5;
  # Another small positive value
  eps2 <- 1e-9;

  # Truncate negative values to a small positive value
  V[V < eps] <- eps;

  # Estimate the initial values of Affinities and Naive Genotyping calls
  WHinit <- WHInit(V[,refs,drop=FALSE]);
  status <- WHinit$status;

  W <- WHinit$W;  # Not really used
  H <- WHinit$H;

  W <- robustWInit(V[,refs,drop=FALSE], H=H);
  H <- robustHInit(V, W=W);

  V <- removeOutliers(V, W=W, H=H);

  # If there is only one allele, no more to do...
  # The algorithm (for one allele) is already a robust estimator
  if (status == 1L || status == 2L) {
    # Shrink average total copy numbers to be close to CN=2.
    totalCNs <- colSums(H[,refs,drop=FALSE]);
    b <- median(totalCNs)/2; # Scale factor
    W <- b*W;
    H <- H/b;
    hasConverged <- NA;
    iter <- NA_integer_;
  } else {
    onesA <- matrix(1, nrow=1L, ncol=I);
    onesB <- matrix(1, nrow=K, ncol=1L);
    ones2 <- matrix(1, nrow=K, ncol=I);

    iter <- 1L;
    hasConverged <- FALSE;
    while (!hasConverged && iter < maxIter) {
      # Remember H from previous iteration to test for convergence
      Hprev <- H;

      # Compute new W solving the system of equations
      H[H < eps] <- eps;
      W <- t(miqr.solve(t(H), t(V)));
      W[W < eps] <- eps;

      # Compute the H
      H <- miqr.solve(W, V);
      H[H < eps] <- eps;

      # Normalizing the W
      norms <- colSums(W);
      norms <- norms + eps2; # Add a small positive value
      W <- W %*% diag(1/norms);
      H <- diag(norms) %*% H;

      # Shrink average total copy numbers to be close to CN=2.
      totalCNs <- colSums(H[,refs,drop=FALSE]);
      b <- median(totalCNs)/2; # Scale factor
      W <- b*W;
      H <- H/b;

      # Converged?
      hasConverged <- (max(abs(Hprev - H)) < acc);

      # Next iteration
      iter <- iter + 1L;
    } # while(...)

    # Robust method for shrinking the average total copy number
    # to close to CN=2.
    Dmat <- rlm(t(H[,refs,drop=FALSE]), matrix(data=2, nrow=ncol(H[,refs,drop=FALSE]), ncol=1L), maxit=maxIterRlm);
    coefs <- Dmat$coefficients;
    H <- diag(coefs) %*% H;
    W <- W %*% diag(1/coefs);

    # Truncate non-positive estimate
    H[H < eps] <- eps;
    W[W < eps] <- eps;
  } # if (status ...)

  # Sanity check (may be removed in the future /HB 2009-03-24)
  stopifnot(nrow(W) == K && ncol(W) == 2L);
  stopifnot(nrow(H) == 2L && ncol(H) == I);

  list(W=W, H=H, hasConverged=hasConverged, nbrOfIterations=iter);
} # fitSnpNmf()


############################################################################
# HISTORY:
# 2010-09-28 [HB]
# o Now argument 'refs' defaults to NULL (not 0), which means all samples.
# o Clean up and robustification of recent edits.
# 2010-06-04 [MO]
# o Added refs as argument.
# 2010-05-18 [MO]
# o Added maxIterRlm as argument.
# 2009-11-18 [HB]
# o Removed internal save() in fitSnpNmf().
# 2009-03-24 [HB]
# o Renamed from Nmf() to fitSnpNmf().  The former name was to generic
#   while our algorithm is rather specific to SNP data.
# o Added optional arguments and internal "constants".
# o Added Rdoc comments.
# o Cleanup.
# 2009-02-15 [MO]
# o Robust method to get the H close to copy number equal to 2.
# 2009-02-05 [MO]
# o Clean the code
# 2009-02-04 [MO]
# o Comment of the lines which try to get the columns of W to be similar
# 2009-01-30 [MO]
# o Robust estimation of W
# o Robust estimation of H
# o With the robust estimations no need to differenciate between status
# o W and H using systems of equations
# o Normalization of the columns of W in each iteration
# o Normalization of the columns of H close to two
############################################################################
