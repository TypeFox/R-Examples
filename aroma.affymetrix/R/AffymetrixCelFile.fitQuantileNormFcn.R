###########################################################################/**
# @set "class=AffymetrixCelFile"
# @RdocMethod fitQuantileNormFcn
#
# @title "Fits quantile normalization functions for the arrays in the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{yTarget}{The target probe signals.}
#   \item{subset}{An optional @numeric @vector specifying the indices of the
#      subset of probes to be used to fit the normalization function.}
#   \item{spar, nknots}{Control parameters passed to
#      @see "stats::smooth.spline".}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a normalization @function.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitQuantileNormFcn", "AffymetrixCelFile", function(this, yTarget, subset=NULL, ..., spar=NULL, nknots=1024, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'yTarget':
  if (!identical(attr(yTarget, "isSorted"), TRUE)) {
    # Sort target signals
    verbose && enter(verbose, "Sorting ", length(yTarget), " target signals");
    yTarget <- sort(yTarget);
    attr(yTarget, "isSorted") <- TRUE;
    verbose && exit(verbose);
  }

  # Argument 'controlParams':
#  spar <- controlParams$spar;
#  nknots <- controlParams$nknots;


  verbose && enter(verbose, "Fitting (quantile) normalization function");

  # Read the probe intensities
  y <- getData(this, fields="intensities", verbose=less(verbose));

  # Sort signals
  verbose && enter(verbose, "Sorting probe signals");
  y <- sort(y);
  verbose && exit(verbose);

  # Fit normalization function
  verbose && enter(verbose, "Fitting smooth spline");
  ok <- !is.na(yTarget) & !is.na(y);
  sp <- smooth.spline(x=y[ok], y=yTarget[ok], spar=spar, nknots=nknots,
                                                          keep.data=FALSE);
  verbose && exit(verbose);

  # Create a minimal 'smooth.spline' object for prediction.
  # Note: You can not do predict(sp$fit, ...) because then there
  # will be a problem with Recall().  See r-devel on 2006-04-05.
  fit <- structure(list(fit=sp$fit), class=class(sp))

  env <- new.env(parent=baseenv());
  assign("fit", fit, envir=env);
  # Create transformation function
  fcn <- function(x, ...) {
    stats::predict(fit, x, ...)$y;
  }
  environment(fcn) <- env;

  verbose && exit(verbose);

  fcn;
}, private=TRUE) # fitQuantileNormFcn()



############################################################################
# HISTORY:
# 2007-06-11
# o BUG FIX: Parameters 'spar' and 'nknots' where retrieved as
#   'controlParams$spar' and 'controlParams$nknots' but 'controlParams'
#   was non-existing.
# 2006-09-14
# o Updated to the new package API.
# 2006-08-25
# o Renamed to normalizeQuantiles().
# o Move to class AffymetrixCelFile and output is now CEL files only.
# 2006-07-21
# o Added more verbose output for normalizeQuantile().
# o BUG FIX: typo for normalizeQuantile(..., format="cel").
# 2006-07-08
# o Added argument 'format' to normalizeQuantile().
# o Added support to write normalized data to CEL files.  Currently we do
#   this by copying the existing CEL file and updating that.
# 2006-05-15
# o Created from AffymetrixDataFile.R.
# 2006-04-05
# o BUG FIX:  fitQuantileNormFcn() returned a function that when called with
#   x values forcing extrapolation, error "Error in Recall(object, xrange) :
#   couldn't find function "predict.smooth.spline.fit" would be thrown.
#   This is because you cannot do predict(sp$fit, ...) but only
#   predict(sp, ...).  Why I don't really know; probably about namespaces.
# 2006-03-18
# o Added argument 'subset' to fitQuantileNormFcn().
# 2006-03-03
# o When creating transformation function in, say, fitQuantileNormFcn(), it
#   is important to create an empty environment for the function otherwise
#   all arguments in the calling function is included too.
# o Added writeApd().  For now it can only write 'intensities'.
############################################################################
