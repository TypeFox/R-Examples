###########################################################################/**
# @set "class=CopyNumberRegions"
# @RdocMethod findAtomicAberrations
# @alias findAtomicAberrations
#
# @title "Finds all possible atomic regions"
#
# \description{
#   @get "title" of a certain length.
# }
#
# @synopsis
#
# \arguments{
#   \item{cnr}{The segments defining the partitioning of the data.}
#   \item{data}{The data  used to test for equality.}
#   \item{H}{A positive @integer specifying how many segments each
#      atomic abberation should contain.}
#   \item{alpha}{A @double in [0,1] specifying the significance level
#      for testing the null-hypothesis that the flanking segments
#      are equal.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @data.frame with K rows, where K >= 0 is the number
#   of atomic aberrations found.
# }
#
# \details{
#   An \emph{aberration of length H} is defined as any H consecutive segments.
#   Each aberrations has two \emph{flanking segments} on each side.
#   Regardless of the content of the aberration, it is possible
#   to test the null-hypothesis that the two flanking segments are
#   equal or not.
#   The two flanking regions are said to be \emph{equal}, if the
#   null-hypothesis of being equal is \emph{not} rejected.
#   If the two flanking regions are called equal, then the contained
#   abberation (of length H) is called \emph{atomic}, otherwise not.
#
#   For consistency one may also define atomic aberrations of length H=0.
#   Consider that an imaginary aberration of zero length splits a single
#   segment into to flanking segments.  Then by construction those two
#   segments are equal.  The case where H=0 is still not implemented.
# }
#
# @examples "../incl/findAtomicAberrations.Rex"
#
# % \references{
# %   [1] \url{http://www.definethat.com/define/7274.htm}
# %   [2] \url{http://en.wikipedia.org/wiki/Atomic_(order_theory)}
# %   [3] \url{http://en.wikipedia.org/wiki/Atomic_(measure_theory)}
# % }
#
# \seealso{
#   ...
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("findAtomicAberrations", "CopyNumberRegions", function(cnr, data, H=1, alpha=0.02, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extractSignals <- function(data, region, ...) {
    data <- extractRegion(data, region=region);
    y <- getSignals(data);
    y <- y[is.finite(y)];
    y;
  } # extractSignals()

  testEquality <- function(dataL, dataR, alpha=0.02, ...) {
    # Test:
    #  H0: muL == muR
    #  H1: muL != muR
    fit <- t.test(dataL, dataR, paired=FALSE, var.equal=TRUE,
                                              alternative="two.sided");
    t <- fit$statistic;
    p <- fit$p.value;
    fit$isSignificant <- (p < alpha);
    isEqual <- (!fit$isSignificant);
    attr(isEqual, "fit") <- fit;

    isEqual;
  } # testEquality()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'data':
  data <- Arguments$getInstanceOf(data, "RawCopyNumbers");

  # Argument 'H':
  H <- Arguments$getInteger(H, range=c(1,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfRegions <- nbrOfRegions(cnr);

  # Nothing to do?
  if (nbrOfRegions < H+2) {
    res <- list(
      atomicRegions=integer(0),
      atomicIslands=integer(0)
    );
    return(res);
  }

  verbose && enter(verbose, "Call equivalent copy-number states by pruning");

  # Initial set of atomic regions
  atomicRegions <- NULL;

  start <- cnr$start;
  stop <- cnr$stop;
  nbrOfPositions <- (nbrOfRegions-H);
  for (rr in 2:nbrOfPositions) {
    verbose && enter(verbose, sprintf("Region #%d of %d", rr, nbrOfPositions));

    # The two flanking regions
    xRangeL <- c(start[rr-1], stop[rr-1]);
    xRangeR <- c(start[rr+H], stop[rr+H]);

    # Extract their data
    dataL <- extractSignals(data, region=xRangeL);
    dataR <- extractSignals(data, region=xRangeR);

    # Test if they are equal
    isEqual <- testEquality(dataL, dataR, alpha=alpha);
    fit <- attr(isEqual, "fit");

    verbose && printf(verbose, "t=%.3f (p=%g), (L==R)=%s\n",
                                       fit$t, fit$p, isEqual);
    # Not needed anymore
    dataL <- dataR <- fit <- NULL;

    # If the two flanking regions are equal, then we have
    # found an atomic region.
    if (isEqual) {
      atomicRegions <- c(atomicRegions, rr);
      verbose && print(verbose, atomicRegions);
    }

    verbose && exit(verbose);
  } # for (rr ...)

  # Table of atomic regions of length K found
  res <- data.frame(
    leftRegion  = atomicRegions-1L,
    rightRegion = atomicRegions+(H-1L)+1L,
    firstRegion = atomicRegions,
    lastRegion  = atomicRegions+(H-1L),
    start       = start[atomicRegions],
    stop        = stop[atomicRegions+(H-1L)]
  );

  # Atomic islands = atomic regions that are not next
  # to another atomic region
  dups <- which(diff(atomicRegions) == 1);
  if (length(dups) > 0) {
    dups <- c(dups, dups+1L);
    atomicIslands <- atomicRegions[-dups];
  } else {
    atomicIslands <- atomicRegions;
  }

  res <- list(
    H=H,
    atomicRegions=atomicRegions,
    atomicIslands=atomicIslands,
    ambigousRegions=setdiff(atomicRegions, atomicIslands),
    res=res
  );

  verbose && exit(verbose);

  res;
}, protected=TRUE) # findAtomicAberrations()


############################################################################
# HISTORY:
# 2010-09-08
# o Added Rdoc comments with an informative code example.
# 2010-09-07
# o Renamed to findAtomicAberrations().
# o Added support for width argument 'H'.
# 2010-07-24
# o CLEAN UP: Now the notation of the code better reflect the algorithm.
# o Now findAtomicRegions() returns ambigous atomic regions too.
# o Added argument 'ylim'.
# 2010-07-20
# o Added argument 'debugPlot'.
# 2010-07-19
# o Added trial version of segmentByPruneCBS().
# o TO DO: Down-weight loci that were close to earlier
#   change points in the succeeding segmentations.
# o Added prototype version of findAtomicRegions().
# o Added prototype version of callByPruning().
# o Created.
############################################################################
