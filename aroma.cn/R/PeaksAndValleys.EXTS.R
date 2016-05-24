###########################################################################/**
# @set "class=PeaksAndValleys"
# @RdocMethod callPeaks
#
# @title "Calls the peaks in peaks-and-valley estimates"
#
# \description{
#   @get "title" to a set of known state.
# }
#
# @synopsis
#
# \arguments{
#  \item{fit}{A KxC @data.frame of peaks-and-valley estimates.}
#  \item{expected}{The expected locations of the peaks to be called.}
#  \item{flavor}{A @character string specifying what flavor of the
#    caller to use.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a Kx(C+2) @data.frame.
# }
#
# \section{Flavors}{
#  If \code{flavor == "all"}, each peak is called to the state with the
#  closest expected value.
#  If \code{flavor == "decreasing"}, the strongest peak is called to the
#  state with the closest expected value, then the second strongest peak
#  is called analogously to one of the remaining states, and so on.
# }
#
# @author "HB, PN"
#
# \seealso{
#   To get peaks-and-valley estimates, use
#   @see "aroma.light::findPeaksAndValleys".
# }
#*/###########################################################################
setMethodS3("callPeaks", "PeaksAndValleys", function(fit, expected=c(-1/2,-1/4,0,+1/4,+1/2)*pi, flavor=c("decreasing", "all"), verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit';
  stopifnot(all(is.element(c("type", "x", "density"), colnames(fit))));

  # Argument 'expected':
  expected <- Arguments$getNumerics(expected);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling peaks");
  verbose && cat(verbose, "Flavor: ", flavor);

  verbose && cat(verbose, "All expected peaks:");
  verbose && print(verbose, expected);

  verbose && enter(verbose, "Extracing peaks");
  subset <- which(fit$type == "peak");
  fitP <- fit[subset,,drop=FALSE];
  verbose && print(verbose, fitP);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calling peaks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  xd <- fitP[,c("x", "density"),drop=FALSE];

  if (flavor == "all") {
    calls <- sapply(xd$x, FUN=function(x) {
      dist <- abs(x - expected);
      which.min(dist);
    });
    r <- seq_along(calls); ## default ranks (used below)
  } else if (flavor == "decreasing") {
    # It is probably better to call the strongest peaks first for which
    # we have more confidence, and then call the other relative to those.
    # /HB 2010-09-19

    # Order peaks by density
    o <- order(xd[,"density"], decreasing=TRUE);
    verbose && cat(verbose, "Reordering:");
    verbose && print(verbose, o);
    # The ranks (for later)
    r <- seq_along(o); r[o] <- r;
    verbose && cat(verbose, "Rank:");
    verbose && print(verbose, r);
    xd <- xd[o,,drop=FALSE];
    verbose && print(verbose, xd);

    # Call the strongest peak first, then the 2nd strongest and so on...
    naValue <- as.integer(NA);
    calls <- rep(naValue, times=nrow(xd));
    expectedLeft <- expected;
    for (kk in seq_len(nrow(xd))) {
      # All expected modes called?
      if (!any(is.finite(expectedLeft))) {
        break;
      }
      # Mode #kk
      x <- xd[kk,"x"];
      dx <- abs(x - expectedLeft);
      call <- which.min(dx);
      expectedLeft[call] <- NA;
      calls[kk] <- call;
    } # for (kk ...)
  } # if (flavor ...)

  verbose && cat(verbose, "Calls:");
  verbose && print(verbose, calls);
  verbose && cat(verbose, "Expected values:");
  verbose && print(verbose, expected[calls]);

  fitC <- cbind(fit, callId=as.integer(NA), call=as.double(NA));
  fitC[subset,"callId"] <- calls[r];
  fitC[subset,"call"] <- expected[calls[r]];
  attr(fitC, "expected") <- expected;

  verbose && print(verbose, fitC);

  verbose && exit(verbose);

  fitC;
}, protected=TRUE) # callPeaks()



##############################################################################
# HISTORY
# 2013-08-04 [HB]
# o CLEANUP: Formally deprecated callPeaks() for data.frame.
# 2012-09-18 [PN]
# o BUG FIX: callPeaks() would return an error when used with flavor "all".
# 2011-10-31 [HB]
# o Added Rdoc comments to callPeaks() for PeaksAndValleys.
# o CLEANUP: Deprecated callPeaks() for data.frame.
# 2010-10-08 [HB]
# o Added callPeaks().
# o Created.
##############################################################################
