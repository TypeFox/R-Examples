###########################################################################/**
# @RdocFunction robustWInit
#
# @title "Robust initialization of the W (affinity) matrix"
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
#  \item{H}{A 2xI @matrix of allele-specific copy-number estimates.}
#  \item{maxIter}{The maximum number of iterations.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a Kx2 @matrix of robustified probe-affinity estimates.
# }
#
# \details{
#   This function utilized a random number generator.
# }
#
# @keyword internal
#*/###########################################################################
robustWInit <- function(V, H, maxIter=50L, ...) {
  # Number of arrays
  I <- ncol(V);
  # Number of probes
  K <- nrow(V);
  # Number of probe pairs
  L <- as.integer(K/2);

  # Sanity check (may be removed in the future /HB 2009-03-24)
  stopifnot(nrow(H) == 2L && ncol(H) == I);

  # A small positive value
  eps <- 1e-5;

  Ws <- matrix(0, nrow=K, ncol=2*maxIter);
  W <- matrix(0, nrow=K, ncol=2L);

  # Create genotyping group of samples
  AA <- which(2*H[2L,] <  H[1L,]);
  BB <- which(H[2L,] > 2*H[1L,]);
  AB <- which(H[2L,] < 2*H[1L,] & 2*H[2L,] >  H[1L,]);

  nAA <- length(AA);
  nAB <- length(AB);
  nBB <- length(BB);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Step 1:
  # In case most of the samples belong to only one group we twist some
  # of them so we "have" signal from both alleles.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  percSamples <- floor(0.85*I);
  if (nAA > percSamples || nAB > percSamples || nBB > percSamples) {
    rrA <- 1:L;
    rrB <- (L+1):K;
    idxs <- sample(I, size=floor(I/2));

    # Majority are heterozygotes?
    if (nAB > percSamples) {
      V[rrA,idxs] <- min(V);
      V[rrB,idxs] <- V[rrB,idxs] * 2;
      H[1L,idxs] <- 0;
      H[2L,idxs] <- 2;
    } else {
      aux <- V[rrA,idxs];
      V[rrA,idxs] <- V[rrB,idxs];
      V[rrB,idxs] <- aux;
      aux <- H[1L,idxs];
      H[1L,idxs] <- H[2L,idxs];
      H[2L,idxs] <- aux;
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Step 2:
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assign arrays into genotype groups (AA, AB, BB).
  AA <- which(2*H[2L,] < H[1L,]);
  BB <- which(H[2L,] > 2*H[1L,]);
  AB <- which(H[2L,] < 2*H[1L,] & 2*H[2L,] >  H[1L,]);

  nAA <- length(AA);
  nAB <- length(AB);
  nBB <- length(BB);

  hasAA <- (nAA > 0L);
  hasAB <- (nAB > 0L);
  hasBB <- (nBB > 0L);

  cont <- 1L;
  for (ii in 1:maxIter) {
    # Select two random samples with different genotypes
    groups <- sample(1:3, size=2L, replace=FALSE);

    sampleAA <- 0L;
    sampleBB <- 0L;
    sampleAB <- 0L;

    # Pick a random AA sample?
    # FIXME: What if hasAA (nAA == 0) is FALSE?!? /HB 2014-04-27
    if(!hasBB || !hasAB || (hasAA && (groups[1L] == 1L || groups[2L] == 1L))) {
      idx <- sample(nAA, size=1L);  # May return 0L? /HB 2014-04-27
      sampleAA <- AA[idx];          # ...which then becomes integer(0)
    }

    # Pick a random BB sample?
    # FIXME: What if hasBB (nBB == 0) is FALSE?!? /HB 2014-04-27
    if(!hasAA || !hasAB || (hasBB && (groups[1L] == 2L || groups[2L] == 2L))) {
      idx <- sample(nBB, size=1L);  # May return 0L? /HB 2014-04-27
      sampleBB <- BB[idx];          # ...which then becomes integer(0)
    }

    # Pick a random AB sample?
    # FIXME: What if hasAB (nAB == 0) is FALSE?!? /HB 2014-04-27
    if(!hasAA || !hasBB || (hasAB && (groups[1L] == 3L || groups[2L] == 3L))) {
      idx <- sample(nAB, size=1L);  # May return 0L? /HB 2014-04-27
      sampleAB <- AB[idx];          # ...which then becomes integer(0)
    }

    # ...and here we get integer(0)*something => integer(0).
    # Comparing (integer(0) > 0L) gives logical(0), which in turn
    # gives an error in the if statements, e.g. if (logical(0)) {}
    # => Error in if (logical(0)) { : argument is of length zero.
    # FIXME: So, the above selection of two random samples is not
    # fully correct/safe.  /HB 2014-04-27
    if (sampleAA*sampleBB > 0L) {
      cc <- c(sampleAA, sampleBB);
    } else if (sampleAB*sampleBB > 0L) {
      cc <- c(sampleAB, sampleBB);
    } else {
      cc <- c(sampleAA, sampleAB);
    }
    dd <- c(cont, cont+1L);
    Ws[,dd] <- t(miqr.solve(t(H[,cc]),t(V[,cc])));

    cont <- cont + 2L;
  } # for (ii ...)

  oddIdxs <- seq(from=1L, to=2L*maxIter, by=2L);
  evenIdxs <- seq(from=2L, to=2L*maxIter, by=2L);
  mediansWA <- rowMedians(Ws[,oddIdxs,drop=FALSE]);
  mediansWB <- rowMedians(Ws[,evenIdxs,drop=FALSE]);

  # Truncate non-positive values
  mediansWA[mediansWA < 0] <- eps;
  mediansWB[mediansWB < 0] <- eps;

  W[,1L] <- mediansWA;
  W[,2L] <- mediansWB;

  # Sanity check (may be removed in the future /HB 2009-03-24)
  stopifnot(nrow(W) == K && ncol(W) == 2L);

  W;
} # robustWInit()


############################################################################
# HISTORY:
# 2014-04-27 [HB]
# o Added comments to highlight potential problems with how the two random
#   samples are choosen in the iteration of the 2nd step.  See also
#   issue report in thread 'CRMA v2 error with Mouse Diversity' on
#   2014-03-25 to the aroma.affymetrix mailing list.  It reports on 'Error
#   in if (sampleAA * sampleBB > 0L) { : argument is of length zero' with
#   "Calls: fit ... fit.ProbeLevelModel -> <Anonymous> -> FUN -> nmfFcn ->
#   robustWInit".
# o CLEANUP: Restructured the iteration of the 2nd step.
# 2009-02-24 [HB]
# o Added Rdoc comments.
# o Cleanig up and standarizing code.
# 2009-02-02 [MO]
# o Change some code to make it more efficient.
# 2009-01-30 [MO]
# o Created
############################################################################
