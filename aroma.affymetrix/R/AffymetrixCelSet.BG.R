###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod calculateParametersGsb
#
# @title "Computes parameters for adjustment of specific binding"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{nbrOfPms}{The number of random PMs to use in estimation.}
#   \item{affinities}{A @numeric @vector of probe affinities.}
#   \item{seed}{An (optional) @integer specifying a temporary random seed
#     to be used during processing.  The random seed is set to its original
#     state when done.  If @NULL, it is not set.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \details{
#   This method is not constant in memory! /HB 2007-03-26
# }
#
# @author "KS"
#*/###########################################################################
setMethodS3("calculateParametersGsb", "AffymetrixCelSet", function(this, nbrOfPms=25000, affinities=NULL, seed=NULL, ..., verbose=FALSE) {
  if (is.null(affinities)) {
    throw("DEPRECATED: Argument 'affinities' to calculateParametersGsb() must not be NULL.");
  }

  args <- list(...);
  for (key in c("path")) {
    if (is.element(key, names(args))) {
      throw(sprintf("DEPRECATED: bgAdjustGcrma() no longer takes argument '%s'.", key));
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  cdf <- getCdf(this);

  verbose && enter(verbose, "Computing parameters for adjustment of specific binding");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
    on.exit(randomSeed("reset"), add=TRUE)
    verbose && printf(verbose, "Random seed temporarily set (seed=%d, kind=\"L'Ecuyer-CMRG\")\n", seed)
  }

  verbose && enter(verbose, "Extracting PM indices");
  cells <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE, verbose=less(verbose,2));
  pmCells <- cells[isPm(cdf, cache=FALSE)];
  # Not needed anymore
  cells <- NULL;
  pmCells <- sort(pmCells);
  verbose && exit(verbose);

  narray <- length(this);

  #  set.seed(1);
  #  was present in original gcrma code; left in here to allow for consistency
  #  check between old and new versions

  # get a sorted random subset of PM to use in parameter estimation
  pmCells.random <- sample(pmCells, size=nbrOfPms);
  pmCells.random <- sort(pmCells.random);
  # Not needed anymore
  pmCells <- NULL;   # Not needed anymore

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "Extracting ", nbrOfPms, " random PM intensities across CEL set");
  # make sure we don't just sample from a single array; avoids problems
  # if we happened to choose a low quality or otherwise aberrant array
  iarray <- sample(seq_len(narray), size=nbrOfPms, replace=TRUE);

  # For each array, read the signals randomized for that array
  # Confirmed to give identical results. /HB 2007-03-26
  pathnames <- getPathnames(this);
  pm.random2 <- vector("double", length=nbrOfPms);
  for (ii in seq_len(narray)) {
    verbose && enter(verbose, sprintf("Array #%d of %d", ii, narray));
    # Cells to be read for this array
    idxs <- which(iarray == ii);
    cells <- pmCells.random[idxs];
    pm.random2[idxs] <- .readCel(pathnames[ii], indices=cells,
                       readIntensities=TRUE, readStdvs=FALSE)$intensities;
    # Not needed anymore
    idxs <- cells <- NULL;
    verbose && exit(verbose);
  } # for (ii ...)
  # Not needed anymore
  iarray <- pathnames <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting probe affinities and fitting linear model")

  aff <- affinities[pmCells.random];
  # Not needed anymore
  pmCells.random <- NULL;  # Not needed anymore

  # Work on the log2 scale
  pm.random2 <- log2(pm.random2);  # Minimize memory usage.
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "Fitting the GCRMA background linear model");
  verbose && str(verbose, pm.random2);
  verbose && str(verbose, aff);
  fit1 <- lm(pm.random2 ~ aff);
  verbose && print(verbose, fit1);
  verbose && exit(verbose);

  verbose && exit(verbose);

  verbose && exit(verbose);

  fit1$coef;
}, private=TRUE)



############################################################################
# HISTORY:
# 2012-11-20 [HB]
# o CLEANUP: Removed bgAdjustRma() for AffymetrixCelSet.
# o CLEANUP: Removed bgAdjustOptical() for AffymetrixCelSet.
# o Now using index 'ii' for arrays everywhere.
# o CLEANUP: bgAdjustOptical(), bgAdjustRma() and bgAdjustGcrma() no
#   longer accepts argument 'path' to be NULL.
# o CLEANUP: bgAdjustGcrma() for AffymetrixCelSet no longer reads
#   and write probe affinity files in the deprecated APD format.
# 2012-10-21 [HB]
# o CLEANUP: Dropped unneeded mkdirs(), because they were all preceeded
#   by an Arguments$getWritablePath().
# 2010-09-29 [HB]
# o ROBUSTNESS: Now bgAdjustGcrma(..., affinities=NULL) is deprecated and
#   throws an exception.
# 2009-08-31 [HB]
# o CLEAN UP: Updated how 'path' is set internally if not specified.
# 2009-04-06 [HB]
# o BUG FIX: The output path of bgAdjustRma() of AffymetrixCelSet would
#   include the full chip type.  It would also give an error if no tags
#   where specified.
# 2009-03-29 [MR]
# o Made slight modifications for bgAdjustGcRma() to work with the
#   newer Gene 1.0 ST arrays.
# 2007-09-06
# o Made calculateParametersGsb() more memory efficient, because it's using
#   the new unlist feature in getCellIndices() of AffymetrixCdfFile.
# 2007-06-30
# o Added .deprecated=TRUE to all methods.
# 2007-03-26
# o Speed up: Using isPm(cdf) instead of readCdfCellIndices() etc.
# o Memory optimization: Found an unlist() without use.names=FALSE.  Saves
#   about 95% of the object, or 250MB of RAM.
# o TO DO: Store GCRMA affinities under annotationData/ and not use the
#   data set paths. /HB
# o Tried to make calculateParametersGsb() constant in memory. Before it
#   read all data for all arrays and then subsampled! /HB
# o BUG FIX: Argument 'minimum' was not used but forced to equal 1.
# 2007-03-22
# o rename gsbParameters to parametersGsb to avoid clash of arguments
#   in bgAdjustGcrma.AffymetrixCelFile().  Not sure why gsbAdjust and
#   gsbParameters are being matched, but there you go.
# 2006-10-10
# o add RMA background correction (normal+exponential)
# 2006-10-06
# o make sure cdf association is inherited
# 2006-10-04
# o Tested, debugged, docs added.
# 2006-09-28
# o Created (based on AffymetrixCelSet.NORM.R).
############################################################################
