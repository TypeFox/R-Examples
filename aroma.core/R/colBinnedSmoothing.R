###########################################################################/**
# @set "class=matrix"
# @RdocMethod colBinnedSmoothing
# @alias colBinnedSmoothing
# @alias binnedSmoothing
# @alias binnedSmoothing.numeric
#
# @title "Binned smoothing of a matrix column by column"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{Y}{A @numeric JxI @matrix (or a @vector of length J.)}
#   \item{x}{A (optional) @numeric @vector specifying the positions of
#     the J entries. The default is to assume uniformly distributed 
#     positions.}
#   \item{w}{A optional @numeric @vector of prior weights for each of 
#     the J entries.}
#   \item{xOut}{Optional @numeric @vector of K bin center locations.}
#   \item{xOutRange}{Optional Kx2 @matrix specifying the boundary
#     locations for K bins, where each row represents a bin \eqn{[x0,x1)}.
#     If not specified, the boundaries are set to be the midpoints
#     of the bin centers, such that the bins have maximum lengths
#     without overlapping.
#     Vice verse, if \code{xOut} is not specified, then \code{xOut} is set
#     to be the mid points of the \code{xOutRange} boundaries.
#   }
#   \item{from, to, by, length.out}{
#     If neither \code{xOut} nor \code{xOutRange} is specified,
#     the \code{xOut} is generated uniformly from these arguments, which 
#     specify the center location of the first and the last bin, and the
#     distance between the center locations, utilizing the
#     @see "base::seq" function.
#     Argument \code{length.out} can be used as an alternative to
#     \code{by}, in case it specifies the total number of bins instead.
#   }
#   \item{FUN}{A @function.}
#   \item{na.rm}{If @TRUE, missing values are excluded, otherwise not.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @numeric KxI @matrix (or a @vector of length K) where
#   K is the total number of bins.
#   The following attributes are also returned:
#   \itemize{
#    \item{\code{xOut}}{The center locations of each bin.}
#    \item{\code{xOutRange}}{The bin boundaries.}
#    \item{\code{count}}{The number of data points within each bin
#         (based solely on argument \code{x}).}
#    \item{\code{binWidth}}{The \emph{average} bin width.}
#  }
# }
#
# \details{
#   Note that all zero-length bins \eqn{[x0,x1)} will get result
#   in an @NA value, because such bins contain no data points.
#   This also means that \code{colBinnedSmoothing(Y, x=x, xOut=xOut)}
#   where \code{xOut} contains duplicated values, will result in
#   some zero-length bins and hence @NA values.
# }
#
# @examples "../incl/colBinnedSmoothing.Rex"
#
# @author
#
# \seealso{
#   @seemethod "colKernelSmoothing".
# }
#
# @keyword array
# @keyword iteration
# @keyword robust
# @keyword univar 
#*/###########################################################################
setMethodS3("colBinnedSmoothing", "matrix", function(Y, x=seq_len(nrow(Y)), w=NULL, xOut=NULL, xOutRange=NULL, from=min(x, na.rm=TRUE), to=max(x, na.rm=TRUE), by=NULL, length.out=length(x), na.rm=TRUE, FUN="median", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'Y'
  n <- nrow(Y);
  k <- ncol(Y);
  
  # Argument 'x'
  if (length(x) != n) {
    throw("Argument 'x' has different number of values than rows in 'Y': ", 
                                                     length(x), " != ", n);
  }

  # Argument 'w'
  if (is.null(w)) {
  } else {
    if (length(w) != n) {
      throw("Argument 'w' has different number of values than rows in 'Y': ", 
                                                       length(w), " != ", n);
    }
  }

  # Argument 'from' & 'to':
  if (is.null(xOut) && is.null(xOutRange)) {
    disallow <- c("NA", "NaN", "Inf");
    from <- Arguments$getNumeric(from, disallow=disallow);
    to <- Arguments$getNumeric(to, range=c(from,Inf), disallow=disallow);
  }

  # Arguments 'by' & 'length.out':
  if (is.null(by) & is.null(length.out)) {
    throw("Either argument 'by' or 'length.out' needs to be given.");
  }
  if (n > 1 && !is.null(by)) {
    by <- Arguments$getNumeric(by, range=c(0,Inf));
  }
  if (!is.null(length.out)) {
    length.out <- Arguments$getInteger(length.out, range=c(0,Inf));
  }

  # Argument 'xOut':
  if (!is.null(xOut)) {
    xOut <- Arguments$getNumerics(xOut);
#    o <- order(xOut);
#    if (!all(diff(o) > 0L)) {
#      throw("Argument 'xOut' must be strictly ordered: ", hpaste(na.omit(xOut)));
#    }
  }

  # Argument 'xOut':
  if (!is.null(xOutRange)) {
    if (!is.matrix(xOutRange)) {
      throw("Argument 'xOutRange' must be a matrix: ", hpaste(class(xOutRange)));
    }
    if (ncol(xOutRange) != 2L) {
      throw("Argument 'xOutRange' must be a matrix with two columns: ", ncol(xOutRange));
    }
    stopifnot(all(xOutRange[,2L] >= xOutRange[,1L]));
  }

  # Arguments 'na.rm':
  na.rm <- Arguments$getLogical(na.rm);

  # Arguments 'FUN':
  if (is.character(FUN)) {
    if (FUN == "median") {
      FUN <- colWeightedMedians;
    } else if (FUN == "mean") {
      FUN <- colWeightedMeans;
    } else {
      throw("Unknown value of argument 'FUN': ", FUN);
    }
  } else if (is.function(FUN)) {
  } else {
    throw("Argument 'FUN' is not a function: ", class(FUN)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Binned smoothing column by column");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Setup (precalculations)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Generate center locations of bins?
  if (is.null(xOut)) {
    # Generate from bin boundaries, or by using seq()?
    if (!is.null(xOutRange)) {
      # Place in the center of the bins
      xOut <- (xOutRange[,1L] + xOutRange[,2L]) / 2;
    } else {
      if (!is.null(by)) {
        xOut <- seq(from=from, to=to, by=by);
      } else {
        xOut <- seq(from=from, to=to, length.out=length.out);
      }
    }
  }

  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);

  # Number of bins
  nOut <- length(xOut);

  # Make sure that 'xOut' (and 'xOutRange') are ordered before
  # binning.  Also make sure to undo afterward.
  o <- order(xOut);
  wasReordered <- (!all(diff(o) > 0L));
  if (wasReordered) {
    ro <- order(o);
    xOut <- xOut[o];
    if (!is.null(xOutRange)) {
      xOutRange <- xOutRange[o,,drop=FALSE];
    }
  }

  # Create 'xOutRange' (or validate existing)
  # [Here 'xOut' must be ordered']
  if (is.null(xOutRange)) {
    # Average bin width
    if (is.null(by)) {
      avgBinWidth <- mean(diff(xOut), na.rm=TRUE);
    } else {
      avgBinWidth <- by;
    }

    # Identify mid points between target locations
    xOutMid <- (xOut[-1L] + xOut[-nOut]) / 2;

    # The width of the first bin is twice the distance between
    # the first and the second 'xOut' such that xOut[1] is the
    # midpoint of (xOutRange[1,1], xOutRange[1,2]).  Likewise
    # for the last bin. /HB 2012-08-26
    xOutMidFirst <- xOutMid[1L] - (xOut[2L] - xOut[1L]);
    xOutMidLast <- xOutMid[nOut-1L] + (xOut[nOut] - xOut[nOut-1L]);
    xOutMid <- c(xOutMidFirst, xOutMid, xOutMidLast);

    naValue <- as.double(NA);
    xOutRange <- matrix(naValue, nrow=nOut, ncol=2);
    xOutRange[,1L] <- xOutMid[-length(xOutMid)];
    xOutRange[,2L] <- xOutMid[-1L];
  } else {
    # Validate
    if (nrow(xOutRange) != nOut) {
      throw("The number of rows in 'xOutRange' does not match the number of bin: ", ncol(xOutRange), " != ", nOut);
    }
  }

  # Assert that the bin boundaries [x0,x1) contains the target bin.
  stopifnot(all(xOutRange[,1L] <= xOut));
  stopifnot(all(xOut <= xOutRange[,2L]));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Smoothing in bins
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Allocate vector of smoothed signals
  naValue <- as.double(NA);
  Ys <- matrix(naValue, nrow=nOut, ncol=k);
  colnames(Ys) <- colnames(Y);

  verbose && enter(verbose, "Estimating signals in each bin");

  verbose && cat(verbose, "Output locations (bin centers):");
  verbose && str(verbose, xOut);

  # Speed up access access
  x0 <- xOutRange[,1L, drop=TRUE];
  x1 <- xOutRange[,2L, drop=TRUE];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Special cases?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  if (identical(FUN, colWeightedMeans) && is.null(w) && require("matrixStats")) {
#    bx <- ...
#    for (cc in seq_len(ncol(Y))) {
#      ys <- matrixStats::binMeans(Y[,cc], x, bx=bx, na.rm=na.rm, count=(cc == 1L));
#      if (cc == 1L) {
#        counts <- attr(ys, "count");
#      }
#      Ys[,cc] <- ys;
#    } # for (cc ...)
#  }

  verbose && cat(verbose, "Summary of bin widths:");
  verbose && print(verbose, summary(x1-x0));

  # Allocate number of counts per bin
  counts <- integer(nOut);

  # For each bin...
  for (kk in seq_len(nOut)) {
    if (kk %% 100 == 0)
      verbose && cat(verbose, kk);

    # Identify data points within the bin
    keep <- which(x0[kk] <= x & x < x1[kk]);
    nKK <- length(keep);
    counts[kk] <- nKK;

    # Nothing to do?
    if (nKK == 0) {
      next;
    }

    # Keep only data points and prior weight within the current bin
    YBin <- Y[keep,,drop=FALSE];
    if (is.null(w)) {
      wBin <- NULL;
    } else {
      wBin <- w[keep];
    }

    value <- FUN(YBin, w=wBin, na.rm=na.rm);

    # Fix: Smoothing over a window with all missing values give zeros, not NA.
    idxs <- which(value == 0);
    if (length(idxs) > 0) {
      # Are these real zeros or missing values?
      YBin <- YBin[idxs,,drop=FALSE];
      YBin <- !is.na(YBin);
      idxsNA <- idxs[colSums(YBin) == 0];    
      value[idxsNA] <- as.double(NA);
    }

#    verbose && str(verbose, value);

    Ys[kk,] <- value;
  } # for (kk ...)

  verbose && exit(verbose);

  # Undo reordering
  if (wasReordered) {
    xOut <- xOut[ro];
    xOutRange <- xOutRange[ro,,drop=FALSE];
    Ys <- Ys[ro,,drop=FALSE];
  }

  # Average bin width
  avgBinWidth <- mean(xOutRange[,2L] - xOutRange[,1L], na.rm=TRUE);

  attr(Ys, "xOut") <- xOut;
  attr(Ys, "xOutRange") <- xOutRange;
  attr(Ys, "counts") <- counts;
  attr(Ys, "binWidth") <- avgBinWidth;

  # Sanity check
  stopifnot(nrow(Ys) == nOut);

  verbose && exit(verbose);

  Ys;
}) # colBinnedSmoothing()



setMethodS3("binnedSmoothing", "numeric", function(y, ...) {
  y <- colBinnedSmoothing(as.matrix(y), ...);
  dim(y) <- NULL;
  y;
})



############################################################################
# HISTORY:
# 2012-08-26
# o BUG FIX: colBinnedSmoothing(..., xOut=xOut) could generate bins
#   at the ends that did not contain the outer most 'xOut' values.
# o BUG FIX: colBinnedSmoothing(..., xOut=xOut) would return binned
#   values in the incorrect order, iff 'xOut' was not ordered.
# o DOCUMENTATION: Clarified that when colBinnedSmoothing() is done over
#   zero-length bins, the output for those bins will be NA.
# 2012-03-14
# o Now colNnnSmoothing() returns a matrix with column name as
#   in argument 'Y'.
# 2012-02-04
# o GENERALIZATION: Now it is possible to call colBinnedSmoothing() with
#   an empty set of input loci, but still requesting a set of output loci,
#   which then will be all missing values.
# o BUG FIX: Argument 'x' of colBinnedSmoothing() would default to the
#   incorrect number of loci.
# 2011-12-15
# o Now colBinnedSmoothing() handles an unordered 'xOut'.
# 2011-12-11
# o ROBUSTNESS: colBinnedSmoothing() now asserts that 'xOut' is ordered.
# 2011-12-10
# o Returning also the bin counts.
# o Now it is possible to fully specify the location and the width
#   of each bin used by colBinnedSmoothing().
# o Moved argument 'xOut' and new 'xOutRange' up front.
# o DOCUMENTATION: Put 'from', 'to', 'by' and 'length.out' into one
#   argument item and reference help for seq() as well.
# o CORRECTNESS: Now bins are by default strictly disjoint, by redefining
#   them as [x0,x1) instead of [x0,x1].
# o CLEANUP: Dropped unused 'xDiff <- xDiff[keep]'.
# 2009-05-16
# o Now colBinnedSmoothing() uses Arguments$getNumerics(), not getDoubles(),
#   where possible.  This will save memory in some cases.
# 2009-05-12
# o Now colBinnedSmoothing() assert that 'from' and 'to' are finite.
# 2009-04-07
# o BUG FIX: When passing a single data points to colBinnedSmoothing(),
#   it would throw the exception: "Range of argument 'by' is out of range
#   [0,0]: [<by>,<by>]".
# 2009-03-23
# o Replace argument 'robust' with more generic 'FUN'.
# 2009-02-11
# o Added more verbose output to colBinnedSmoothing().
# 2009-02-07
# o Created.
############################################################################
