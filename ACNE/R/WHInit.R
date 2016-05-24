###########################################################################/**
# @RdocFunction WHInit
#
# @title "Initialization of the W and H matrices"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{V}{An KxI @matrix where I is the number of arrays and
#    K is the number of probes where K should be even (K=2L).}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @list:
#  \item{W}{A Kx2 @matrix of initial probe-affinity estimates.}
#  \item{H}{A 2xI @matrix of initial allele-specific copy-number estimates.}
#  \item{status}{An @integer specifying the status of the initialization:
#    0=normal case, 1=only one allele (either AA or BB), or
#    2=all samples are AB.
#  }
# }
#
# \details{
#   The allele-specific copy number estimates are estimated using a
#   naive genotyping algorithm.
#   The probe-affinities are estimated using a pseudo inverse.
# }
#
# @keyword internal
#*/###########################################################################
WHInit <- function(V, ...) {
  # Number of arrays
  I <- ncol(V);
  # Number of probes
  K <- nrow(V);
  # Number of probe pairs
  L <- as.integer(K/2);

  # A small positive value
  eps <- 0.0001;

  H <- matrix(0, nrow=2L, ncol=I);
  W <- matrix(0, nrow=K, ncol=2L);
  rrA <- 1:L;
  rrB <- (L+1):K;
  rrBA <- c(rrB, rrA);
  PMA <- V[rrA,,drop=FALSE];
  PMB <- V[rrB,,drop=FALSE];

  # We distinguish three locations:
  #  (1) AA (PMA > 2 PMB),
  #  (2) AB (PMA < 2PMB & PMB < 2PMA), and
  #  (3) BB (PMB > 2PMB).
  # We apply this test for each of the probes and we use majority voting.
  H[1,] <- as.integer(colMeans(PMA > 0.5*PMB) > 0.5);
  H[2,] <- as.integer(colMeans(PMB > 0.5*PMA) > 0.5);

  summary <- 2*H[1L,] + H[2L,];
  dummy <- unique(summary);

  status <- 0L;
  # If all the samples have the same genotype, it is a special case
  if (length(dummy) == 1L) {
    # We have only one Genotype
    # Case AA or BB
    if (prod(H[,1L]) == 0) {
      #print('only one allele AA or BB');
      # Use the median for the first column of W
      W[,1L] <- rowMedians(V)/2;
      # Flip it for the second column
      W[,2L] <- W[rrBA,1L];
      # Change both of them if it was BB
      H <- H*2;
      if (H[2L,1L] == 1){
        # Was it BB?
        W <- W[,c(2L,1L),drop=FALSE];
        H <- H[c(2L,1L),,drop=FALSE];
      }
      status <- 1L;
    } else {
      #disp('only samples AB')
      W[,1L] <- rowMedians(V);
      W[,2L] <- W[,1L];
      # In this case there is no way to compute the cross hybridization
      # We assume that there is no cross hybridization (just to asssume
      # something :-)
      W[rrB,1L] <- eps;
      W[rrA,2L] <- eps;
      status <- 2L;
    }
  } else {
    # Normal case
    aux <- colSums(H);
    aux <- rep(aux, times=2L);
    dim(aux) <- c((length(aux)/2), 2);
    aux <- t(aux);
    H <- 2 * H/aux;
    H[is.na(H)] <- 0;
    W <- t(miqr.solve(t(H),t(V)));
    W[W < 0] <- eps;

    # Sometimes, there are errors in the genotyping... Check correlation
    corDiff <- cor(W[,1L],W[rrBA,2L]) - cor(W[,1L],W[,2L]);
    if (is.na(corDiff) || corDiff < 0.1) {
      #print('Too large Correlation')
      #print('Solving for one allele')
      W0 <- W;
      W[,1L] <- rowMedians(W0);
      W[,2L] <- W0[rrBA,1L];
      H <- miqr.solve(W, V);
      H[H < 0] <- 0;
      status <- 1L;
    }
  }

  # Sanity check (may be removed in the future /HB 2009-03-24)
  stopifnot(nrow(W) == K && ncol(W) == 2L);
  stopifnot(nrow(H) == 2L && ncol(H) == I);

  list(W=W, H=H, status=status);
} # WHInit()

############################################################################
# HISTORY:
# 2009-03-24 [HB]
# o Added Rdoc comments.
# o Cleaned up code.
# 2009-02-02 [MO]
# o Change some code and make it more efficient
# 2009-01-28 [HB]
# o BUG FIX: The code of WHInit() assumed 20 probes in one place.
############################################################################
