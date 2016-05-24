###########################################################################/**
# @set "class=numeric"
# @RdocMethod callXXorXY
# @alias callXXorXY
#
# @title "Calls XX or XY from ChrX allele B fractions of a normal sample"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{betaX}{A @numeric @vector containing ChrX allele B fractions.}
#  \item{betaY}{A optional @numeric @vector containing ChrY allele B fractions.}
#  \item{flavor}{A @character string specifying the type of algorithm used.}
#  \item{adjust}{A postive @double specifying the amount smoothing for
#    the empirical density estimator.}
#  \item{...}{Additional arguments passed to
#    @see "aroma.light::findPeaksAndValleys".}
#  \item{censorAt}{A @double @vector of length two specifying the range
#    for which values are considered finite.  Values below (above) this
#    range are treated as -@Inf (+@Inf).}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a ...
# }
#
# \section{Missing and non-finite values}{
#   Missing and non-finite values are dropped before trying to call XX or XY.
# }
#
# @author "HB, PN"
#
# \seealso{
#   Internally @see "aroma.light::findPeaksAndValleys" is used to identify
#   the thresholds.
# }
#*/###########################################################################
setMethodS3("callXXorXY", "numeric", function(betaX, betaY=NULL, flavor=c("density"), adjust=1.5, ..., censorAt=c(-0.5,+1.5), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'betaX':
  betaX <- as.double(betaX);

  # Argument 'betaY':
  if (!is.null(betaY)) {
    betaY <- as.double(betaY);
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'adjust':
  adjust <- as.double(adjust);
  if (length(adjust) != 1) {
    stop("Argument 'adjust' must be single value: ", adjust);
  }
  if (adjust <= 0) {
    stop("Argument 'adjust' must be positive: ", adjust);
  }

  # Argument 'censorAt':
  censorAt <- Arguments$getDoubles(censorAt, length=c(2,2));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling gender from allele B fractions (BAFs)");

  betaT <- betaX;
  betaT[betaT < censorAt[1]] <- -Inf;
  betaT[betaT > censorAt[2]] <- +Inf;
  betaT <- betaT[is.finite(betaT)];
  fit <- .findPeaksAndValleys(betaT, adjust=adjust, ...);
  isXYByChrX <- (sum(fit$type == "peak") == 2);

  if (isXYByChrX && (length(betaY) > 100)) {
    betaT <- betaY;
    betaT[betaT < censorAt[1]] <- -Inf;
    betaT[betaT > censorAt[2]] <- +Inf;
    fit <- .findPeaksAndValleys(betaT, adjust=adjust, ...);
    isXYByChrY <- (sum(fit$type == "peak") == 2);
    if (isXYByChrY != isXYByChrX) {
      throw("Allele B fractions for ChrX and ChrY are inconsistent.");
    }
  }

  res <- ifelse(isXYByChrX, "XY", "XX");

  verbose && exit(verbose);

  res;
}) # callXXorXY()


###########################################################################
# HISTORY:
# 2012-04-16 [HB]
# o Now callXXorXY() explicitly requires aroma.light.
# 2011-03-03 [HB]
# o TYPO: Used betaX[is.finite(betaT)] instead of betaT[is.finite(betaT)],
#   but the results would have been identical either way.
# 2010-07-22 [PN]
# o No longer calling gender from chr Y when gender is estimated as
# "XX" from chr X.
# 2009-11-03
# o Created.
###########################################################################
