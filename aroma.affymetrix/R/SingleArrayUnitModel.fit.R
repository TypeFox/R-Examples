###########################################################################/**
# @set "class=SingleArrayUnitModel"
# @RdocMethod fit
#
# @title "Estimates the model parameters"
#
# \description{
#  @get "title" for all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{The arrays to be fitted.
#     If @NULL, all arrays are considered.
#     If \code{remaining}, only non-fitted arrays are considered.
#   }
#   \item{units}{The units to be fitted.
#     If @NULL, all units are considered.
#     If \code{remaining}, only non-fitted units are considered.
#   }
#   \item{...}{Arguments passed to @seemethod "readUnits".}
#   \item{force}{If @TRUE, already fitted units are re-fitted, and
#     cached data is re-read.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns nothing.
# }
#
# \details{
#  All estimates are stored to file.
#
#  The parameter estimates specific to each array,
#  typically "chip effects",
#  are stored in array specific files.
#
#   Array-specific estimates [K = nbr of arrays]:
#    theta [K doubles] (chip effects), sd(theta) [K doubles],
#    isOutlier(theta) [K logicals]
#
#   For each array and each unit group, we store:
#     1 theta, 1 sd(theta), 1 isOutlier(theta), i.e. (float, float, bit)
#   => For each array and each unit (with \eqn{G_j} groups), we store:
#     \eqn{G_j} theta, \eqn{G_j} sd(theta), \eqn{G_j} isOutlier(theta),
#   i.e. \eqn{G_j}*(float, float, bit).
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fitOneArray", "SingleArrayUnitModel", function(this, array="remaining", units="remaining", ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the some basic information about this model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  nbrOfArrays <- length(ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  array <- Arguments$getIndex(array, max=nbrOfArrays);

  # Argument 'units':
  doRemaining <- FALSE;
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  } else if (identical(units, "remaining")) {
    doRemaining <- TRUE;
  } else {
    throw("Unknown mode of argument 'units': ", mode(units));
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model-fit function
  fitUnit <- getFitUnitFunction(this);

  ces <- getChipEffectSet(this);

  df <- ds[[array]];
  cef <- ces[[array]];

  # Sanity check
  stopifnot(identical(getFullName(df), getFullName(cef)));

  name <- getName(df);
  verbose && enter(verbose, sprintf("Fitting array #%d ('%s')", array, name));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify units to be fitted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    units <- 1:nbrOfUnits;
  } else if (doRemaining) {
    verbose && enter(verbose, "Identifying non-estimated units")
    units <- findUnitsTodo(cef, verbose=less(verbose));
    nbrOfUnits <- length(units);
    verbose && exit(verbose);
  } else {
    # Fit only unique units
    units <- unique(units);
    nbrOfUnits <- length(units);
  }
  verbose && printf(verbose, "Getting model fit for %d units.\n", nbrOfUnits);

  # Identify which of the requested units have *not* already been estimated
  if (!doRemaining) {
    if (force) {
      verbose && printf(verbose, "All of these are forced to be fitted.\n");
    } else {
      units <- findUnitsTodo(cef, units=units, verbose=less(verbose));
      nbrOfUnits <- length(units);
      verbose && printf(verbose, "Out of these, %d units need to be fitted.\n", nbrOfUnits);
    }
  }

  # Nothing more to do?
  if (nbrOfUnits == 0) {
    verbose && cat(verbose, "Nothing to do.");
    verbose && exit(verbose);
    return(NULL);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Timers BEGIN
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  startTime <- processTime();
  timers <- list(total=0, read=0, fit=0, writePaf=0, writeCes=0, gc=0);
  tTotal <- processTime();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the CEL intensities by units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Units: ");
  verbose && str(verbose, units);

  tRead <- processTime();
  y <- readUnits(this, array=array, units=units, ..., force=force,
                                     cache=FALSE, verbose=less(verbose));
  timers$read <- timers$read + (processTime() - tRead);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting model");
  tFit <- processTime();

  fit <- lapply(y, FUN=fitUnit);

  timers$fit <- timers$fit + (processTime() - tFit);
  y <- NULL; # Not needed anymore (to minimize memory usage)
  verbose && str(verbose, fit[1]);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store chip-effect estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Storing chip-effect estimates");
  tWriteCes <- processTime();
  updateUnits(cef, units=units, data=fit, verbose=less(verbose));
  timers$writeCes <- timers$writeCes + (processTime() - tWriteCes);
  verbose && exit(verbose);

  fit <- NULL; # Not needed anymore

  # Garbage collection
  tGc <- processTime();
  gc <- gc();
  timers$gc <- timers$gc + (processTime() - tGc);
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Timers END
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  timers$total <- timers$total + (processTime() - tTotal);

  totalTime <- processTime() - startTime;
  if (verbose) {
    nunits <- length(units);
    # Get distribution of what is spend where
    timers$write <- timers$writePaf + timers$writeCes;
    t <- lapply(timers, FUN=function(timer) unname(timer[3]));
    t <- unlist(t);
    t <- 100 * t / t["total"];
    printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%% (of which %.2f%% is for encoding/writing chip-effects), Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], 100*t["writeCes"]/t["write"], t["gc"]);
  }

  verbose && exit(verbose);

  invisible(units);
}, protected=TRUE);  # fitOneArray()



setMethodS3("fit", "SingleArrayUnitModel", function(this, arrays=NULL, units="remaining", ..., force=FALSE, verbose=FALSE) {


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit model array by array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (aa in seq_along(arrays)) {
    array <- arrays[aa];
    verbose && enter(verbose, aa);

    fitOneArray(this, array=array, units=units, verbose=verbose);

    verbose && exit(verbose);
  } # for (aa ...)


  invisible(units);
})




############################################################################
# HISTORY:
# 2008-09-03
# o Created from ProbeLevelModel.fit.R.
############################################################################
