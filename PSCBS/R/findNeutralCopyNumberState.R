###########################################################################/**
# @RdocDefault findNeutralCopyNumberState
#
# @title "Call segments to be copy neutral based on allelic imbalance calls and total copy number estimates"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{C}{A @numeric @vector of region-level total copy number estimates.}
#   \item{isAI}{A @logical @vector of "allelic imbalance" calls.}
#   \item{weights}{An optional @numeric @vector of non-negative weights.}
#   \item{...}{Further argumants to be passed to the density estimation
#     function.}
#   \item{minDensity}{A @numeric value, below which density peaks are
#     discarded.}
#   \item{flavor}{A @character string specifying how to identify the
#     mode of the AB segments.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#   A @logical @vector of "neutral copy number state" calls.
# }
#
# @author "PN, HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("findNeutralCopyNumberState", "default", function(C, isAI, weights=NULL, ..., minDensity=1e-10, flavor=c("firstPeak", "maxPeak"), verbose=FALSE) {
  # This will load the 'aroma.light' namespace, if not already done.
  findPeaksAndValleys <- .use("findPeaksAndValleys", package="aroma.light");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'C':
  C <- Arguments$getNumerics(C);
  nbrOfLoci <- length(C);

  # Argument 'isAI':
  length2 <- rep(nbrOfLoci, times=2);
  isAI <- Arguments$getLogicals(isAI, length=length2, disallow=NULL);

  # Argument 'weights':
  if (!is.null(weights)) {
    weights <- Arguments$getNumerics(weights, range=c(0, Inf), length=length2);
  }

  # Argument 'minDensity':
  minDensity <- Arguments$getDouble(minDensity);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Identifying segments that are copy neutral states");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify segments in allelic balance
  isAB <- !isAI;

  # Identify segments that cannot be called
  isNA <- (is.na(isAB) | is.na(C));

  # Only segments in allelic balance can be considered to be neutral
  isNeutral <- isAB;

  # Extracting segments in allelic balance
  idxs <- which(isAB);
  n <- length(idxs);
  verbose && cat(verbose, "Number of segments in allelic balance: ", n);

  # Special cases?
  if (n == 0) {
    # No segments are in allelic balance
    verbose && exit(verbose);
    return(isNeutral);
  } else if (n == 1) {
    # Only one segment is in allelic balance.  The best we can do
    # is to call that segment neutral.
    verbose && exit(verbose);
    return(isNeutral);
  } else if (n < 5) {
    # What to do when the number of segments is really low? /HB 2010-09-09
    warning("The calling of regions in a copy-neutral state is uncertain, because there are less than five (5) regions in allelic balance: ", n);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Look only segments in allelic balance
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset and standardize weights
  if (!is.null(weights)) {
    weights <- weights[idxs];
    weights <- weights / sum(weights);
  }
  y <- C[idxs];
  idxs <- NULL; # Not needed anymore

  if (verbose) {
    cat(verbose, "Data points:");
    df <- data.frame(C=y, weights=weights);
    print(verbose, head(df));
    str(verbose, df);
    df <- NULL; # Not needed anymore
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate the empirical density
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit <- findPeaksAndValleys(y, weights=weights, ...);
  verbose && cat(verbose, "Fit:");

  verbose && cat(verbose, "Fit filtered by 'minDensity':");
  ok <- (fit[,"density"] > minDensity);
  verbose && print(verbose, fit[ok,]);

  # Look for peaks with enough density
  isPeak <- (fit[,"type"] == "peak") & ok;
  idxs <- which(isPeak);

  # Sanity check
  stopifnot(length(idxs) >= 1);

  # Extract the first peak
  if (flavor == "firstPeak") {
    idx <- idxs[1];
  } else if (flavor == "maxPeak") {
    idx <- idxs[which.max(fit[idxs,"density"])];
  }

  neutralC <- fit[idx,"x"];

  verbose && cat(verbose, "Neutral copy number:");
  verbose && cat(verbose, "Mode at: ", neutralC);
  verbose && cat(verbose, "Mode ampliture: ", fit[idx,"density"]);

  # If there is more than one peak, we should only call segments that
  # are not part of that other peak.
  if (idx+1 <= nrow(fit)) {
    nextValleyC <- fit[idx+1, "x"];
  } else {
    nextValleyC <- Inf;
  }
  verbose && cat(verbose, "Upper range at: ", nextValleyC);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call copy-neutral regions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isNeutral <- isNeutral & (C < nextValleyC);

  # Segments with missing values cannot be called
  isNeutral[isNA] <- NA;

  verbose && cat(verbose, "Neutral region calls:");
  verbose && summary(verbose, isNeutral);

  verbose && exit(verbose);

  isNeutral;
}) # findNeutralCopyNumberState()


##############################################################################
# HISTORY
# 2013-09-26 [HB]
# o CLEANUP: Now findNeutralCopyNumberState() no longer attached
#   'aroma.light', but only loads its namespace.
# 2013-03-19 [HB]
# o Added argument 'flavor' to findNeutralCopyNumberState() specifying how
#   to identify the main mode of the AB segments.
# 2012-02-24 [HB]
# o Moved findNeutralCopyNumberState() from aroma.light.
# 2012-02-23 [HB]
# o Renamed argument 'densityThreshold' to 'minDensity'.
# 2011-07-10 [HB]
# o Made findNeutralCopyNumberState() a default method.
# o Made the Rd help "internal".
# 2010-09-09 [HB]
# o Now segments with missing values are not called.
# o Added support for the case when there is no peak/no segments in AB.
# o Added support for the case when there is only one weak.
# o Added sanity checks.
# 2010-09-08 [PN]
# o Created.
##############################################################################
