###########################################################################/**
# @RdocClass DChipQuantileNormalization
#
# @title "The DChipQuantileNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a special @see "QuantileNormalization"
#  using smooth-splines.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#       @see "QuantileNormalization".}
#   \item{robust}{If @TRUE, the normalization function is estimated
#       robustly, otherwise not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   This normalization method implements the two-pass algorithm described
#   in Bengtsson et al. (2008).
# }
#
# \references{
#   [1] H. Bengtsson, R. Irizarry, B. Carvalho, & T.P. Speed.
#       Estimation and assessment of raw copy numbers at the single
#       locus level, Bioinformatics, 2008.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("DChipQuantileNormalization", function(..., robust=FALSE) {
  # Arguments 'robust':
  robust <- Arguments$getLogical(robust);

  extend(QuantileNormalization(...), "DChipQuantileNormalization",
    .robust = robust,
    .exclCells = NULL
  );
})


setMethodS3("as.character", "DChipQuantileNormalization", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  nExcl <- length(getExclCells(this));
  n <- nbrOfCells(getCdf(getInputDataSet(this)));
  s <- c(s, sprintf("Number of cells excluded (when fitting): %d (%.1f%%)",
                                                         nExcl, 100*nExcl/n));
  s;
}, protected=TRUE)


setMethodS3("getParameters", "DChipQuantileNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  params$robust <- this$.robust;
  subsetToAvg <- params$subsetToAvg;

  exclCells <- getExclCells(this);
  if (!is.null(exclCells)) {
    if (is.null(subsetToAvg)) {
      ds <- getInputDataSet(this);
      cdf <- getCdf(ds);
      subsetToAvg <- seq_len(nbrOfCells(cdf));
      subsetToAvg <- setdiff(subsetToAvg, exclCells);
      subsetToAvg <- sort(subsetToAvg);
    }
  }
  params$subsetToAvg <- subsetToAvg;

  params;
}, protected=TRUE)


setMethodS3("getSubsetToUpdate", "DChipQuantileNormalization", function(this, ...) {
  subsetToUpdate <- NextMethod("getSubsetToUpdate");
  if (is.null(subsetToUpdate)) {
    if (is.null(this$.typesToUpdate)) {
    } else if (this$.typesToUpdate == "pm") {
    }
    this$.subsetToUpdate <- subsetToUpdate;
  }
  subsetToUpdate;
}, private=TRUE)


setMethodS3("getExclCells", "DChipQuantileNormalization", function(this, ..., verbose=FALSE) {
  this$.exclCells;
}, protected=TRUE)


setMethodS3("addExclCells", "DChipQuantileNormalization", function(this, cells, ..., verbose=FALSE) {
  cells <- c(this$.exclCells, cells);
  cells <- unique(cells);
  cells <- sort(cells);
  this$.exclCells <- cells;
  invisible(this);
}, protected=TRUE)


setMethodS3("excludeChrXFromFit", "DChipQuantileNormalization", function(this, ..., verbose=FALSE) {
  # Identify cells on chromosome X
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  gi <- getGenomeInformation(cdf);
  units <- getUnitsOnChromosome(gi, 23);
  cells <- getCellIndices(cdf, units=units, useNames=FALSE, unlist=TRUE);

  # Add them to the list of cells to be excluded
  addExclCells(this, cells);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to
#       @see "aroma.light::normalizeQuantileSpline.numeric".}
#   \item{force}{If @TRUE, data already normalized is re-normalized,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "DChipQuantileNormalization", function(this, ..., force=FALSE, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Quantile normalizing (using smooth splines) data set");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  if (!force && isDone(this)) {
#    verbose && cat(verbose, "Already normalized");
#    verbose && exit(verbose);
#    outputDataSet <- getOutputDataSet(this);
#    return(invisible(outputDataSet));
#  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit non-robustly (faster and more memory efficient).
  robust <- this$.robust;

  # Get input data set
  ds <- getInputDataSet(this);

  cdf <- getCdf(ds);

  # Get (and create) the output path
  outputPath <- getPath(this);

  params <- getParameters(this);

  # Retrieve/calculate the target distribution
  xTarget <- getTargetDistribution(this, verbose=less(verbose));
  xTarget <- sort(xTarget, na.last=TRUE);
  verbose && cat(verbose, "Target distribution: ");
  verbose && str(verbose, xTarget);

  # TO DO: Add function to "expand" 'xTarget' if of different length
  # than 'x' and 'w'.  /HB 2007-04-11. DONE 2008-02-23.
  if (length(xTarget) != nbrOfCells(cdf)) {
    # See normalizeQuantileSpline() for why this is ok/the way to do it.
    xTarget <- c(xTarget, rep(NA, nbrOfCells(cdf)-length(xTarget)));
    verbose && cat(verbose, "Expanded target distribution (now with NAs): ");
    verbose && str(verbose, xTarget);
  }

  # Get algorithm parameters
  verbose && cat(verbose, "typesToUpdate: ");
  verbose && str(verbose, params$typesToUpdate);
  verbose && cat(verbose, "subsetToUpdate: ");
  verbose && str(verbose, params$subsetToUpdate);
  subsetToUpdate <- identifyCells(cdf, indices=params$subsetToUpdate,
                         types=params$typesToUpdate, verbose=less(verbose));
  verbose && str(verbose, subsetToUpdate);

  # Exclude certain cells when *fitting* the normalization function
  excl <- getExclCells(this);
  if (length(excl) > 0) {
    verbose && enter(verbose, "Excluded some cells when fitting normalization function");

    w <- rep(1, nbrOfCells(cdf));
    w[excl] <- 0;

    # If not all cells, get weights in the same order as the data points 'x'.
    if (!is.null(subsetToUpdate)) {
      w <- w[subsetToUpdate];
    }

    # Standardize weights to sum to one.
    w <- w / sum(w, na.rm=TRUE);
    verbose && printf(verbose, "Cell weights (sum = %.2f):\n",
                                                     sum(w, na.rm=TRUE));
    verbose && summary(verbose, w);
    verbose && exit(verbose);
  } else {
    w <- NULL;
  }


  # Not needed anymore
  # Not needed anymore
  excl <- NULL;

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing ", length(ds), " arrays");
  dataFiles <- list();
  for (kk in seq_along(ds)) {
    verbose && enter(verbose, "Array #", kk);
    df <- ds[[kk]];
    verbose && print(verbose, df);

    filename <- basename(getPathname(df));
    filename <- gsub("[.]cel$", ".CEL", filename);  # Only output upper case!
    pathname <- Arguments$getWritablePathname(filename, path=outputPath);
    pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);

    # Already normalized?
    if (skip && isFile(pathname)) {
      verbose && cat(verbose, "Normalized data file already exists: ",
                                                                   pathname);
      # CDF inheritance
      dataFiles[[kk]] <- fromFile(df, pathname);
      verbose && exit(verbose);
      next;
    }

    # Get all probe signals
    verbose && enter(verbose, "Reading probe intensities");
    x <- getData(df, indices=subsetToUpdate, fields="intensities",
                                        verbose=less(verbose,2))$intensities;
    verbose && str(verbose, x);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # TO DO: Add function to "expand" 'xTarget' if of different length
    # than 'x' and 'w'.  /HB 2007-04-11, 2008-02-23

    x <- .normalizeQuantileSpline(x, w=w, xTarget=xTarget,
                                       sortTarget=FALSE, robust=robust, ...);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Write normalized data to file
    verbose && enter(verbose, "Writing normalized probe signals");

    # Write to a temporary file (allow rename of existing one if forced)
    isFile <- (!skip && isFile(pathname));
    pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing");
    createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose));
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing normalized intensities");
    .updateCel(pathnameT, indices=subsetToUpdate, intensities=x);

    # Not needed anymore
    x <- NULL;

    # Rename temporary file
    popTemporaryFile(pathnameT, verbose=verbose);

    verbose && exit(verbose);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Return new normalized data file object
    dataFiles[[kk]] <- fromFile(df, pathname);

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Not needed anymore
  # Not needed anymore
  w <- NULL;

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  # Create result set
  outputDataSet <- newInstance(ds, dataFiles);
  setCdf(outputDataSet, cdf);

  # Update the output data set
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);

  outputDataSet;
})


############################################################################
# HISTORY:
# 2008-02-23
# o BUG FIX: When excluding cells from the fit, we would get an error saying
#   the length of the target distribution is not the same as the data to
#   be normalized. I had put this up on the todo list already 2007-04-11,
#   but it is first now I got around to fix it.
# 2007-09-06
# o Made excludeChrXFromFit() more memory efficient, because it's using
#   the new unlist feature in getCellIndices() of AffymetrixCdfFile.
# 2007-04-08
# o Added argument 'robust' to the constructor.
# 2007-03-28
# o Added getParameters() so that excluded cells are also excluded when
#   the target distribution is calculated by the average.
# 2007-03-22
# o Added test code for excluding some cells by giving them weight zero
#   when fitting the normalization function.
# o TO DO: Add code to exclude non-diploid data points from the estimation
#   of the normalization function, e.g. male chromosome X signals should
#   not included if the target distribution was calculated for females only.
#   To simplify it, we could exclude all chromosome X signals regardless
#   of ploidy.  In that we don't have to know the ploidy. The those should
#   also be excluded when estimating the target distribution.
#   Potential problems: If the majority of the signals are from chrX, then
#   it does not work.
# 2006-12-11
# o Created to immitate the oligo package.
############################################################################
