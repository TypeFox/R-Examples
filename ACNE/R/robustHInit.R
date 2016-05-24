###########################################################################/**
# @RdocFunction robustHInit
#
# @title "Robust initialization of the H (copy number) matrix"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{V}{An KxI @matrix where I is the number of arrays and K is the
#     number of probes where K should be even (K=2L).}
#  \item{W}{A Kx2 @matrix of probe-affinity estimates.}
#  \item{maxIter}{The maximum number of iteration.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a 2xI @matrix of robustified allele-specific copy-number estimates.
# }
#
# \details{
#   This function utilized a random number generator.
# }
#
# @keyword internal
#*/###########################################################################
robustHInit <- function(V, W, maxIter=5L, ...) {
  # Number of arrays
  I <- ncol(V);
  # Number of probes
  K <- nrow(V);

  # Sanity check (may be removed in the future /HB 2009-03-24)
  stopifnot(nrow(W) == K && ncol(W) == 2L);

  # V = W * H

  #[Vi1... Vin; Vj1... Vjn] = W * H
  # 2 random probes a especific number of iterations

  H <- matrix(0, nrow=2L, ncol=I);
  Haux <- matrix(0, nrow=maxIter*K, ncol=I);

  contHaux <- 1L;
  for (ii in 1:maxIter) {
    probes <- sample(K);
    oddIdxs <- seq(from=1L, to=length(probes/2), by=2L);
    for (jj in oddIdxs) {
      pp <- c(probes[jj], probes[jj+1]);
      Haux[c(contHaux,contHaux+1L),] <- miqr.solve(W[pp,],V[pp,]);
      contHaux <- contHaux + 2L;
    } # for (jj ...)
  } # for (ii ...)

  # Truncate non-positive values
  Haux[Haux < 0] <- 0;

  oddIdxs <- seq(from=1L, to=maxIter*K, by=2L);
  evenIdxs <- seq(from=2L, to=maxIter*K, by=2L);
  H[1L,] <- colMedians(Haux[oddIdxs,,drop=FALSE]);
  H[2L,] <- colMedians(Haux[evenIdxs,,drop=FALSE]);

  # Sanity check (may be removed in the future /HB 2009-03-24)
  stopifnot(nrow(H) == 2L && ncol(H) == I);

  H;
} # robustHInit()


############################################################################
# HISTORY:
# 2009-03-24 [HB]
# o Renamed from RobustHinit() to robustHinit().
# o Cleaned up code.
# o Added Rdoc comments.
# 2009-02-02 [MO]
# o Change some code to make more efficient and change the name of the
#   indexes.
# 2009-01-30 [MO]
# o Created
############################################################################
