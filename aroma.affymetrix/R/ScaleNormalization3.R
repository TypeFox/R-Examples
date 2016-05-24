###########################################################################/**
# @RdocClass ScaleNormalization3
#
# @title "The ScaleNormalization3 class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization function that transforms the
#  probe-level signals towards the same scale.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform3".}
#   \item{targetAvg}{A @numeric value.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("ScaleNormalization3", function(..., targetAvg=4400) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  targetAvg <- Arguments$getDouble(targetAvg, range=c(1,Inf));

  extend(ProbeLevelTransform3(...), "ScaleNormalization3",
    .targetAvg = targetAvg
  )
})



setMethodS3("getParameters", "ScaleNormalization3", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    targetAvg = this$.targetAvg
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)



setMethodS3("fitOne", "ScaleNormalization3", function(this, df, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  df <- Arguments$getInstanceOf(df, "AffymetrixCelFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Fitting normalization function for one array");
  verbose && cat(verbose, "Full name: ", getFullName(df));

  verbose && enter(verbose, "Getting algorithm parameters");
  params <- getParameters(this, expand=TRUE);
  cells <- params$cellsToFit;
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  shift <- params$shift;
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading signals");
#  y <- extractMatrix(df, units=cells, drop=TRUE, verbose=verbose);  # NOTE: 'units' :(
  y <- getData(df, indices=cells, fields="intensities", drop=TRUE, verbose=verbose);
  verbose && str(verbose, y);
  verbose && exit(verbose);

  # Shift?
  if (shift != 0) {
    verbose && enter(verbose, "Shifting signals");
    y <- y + shift;
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Estimating mean parameter");
  mu <- median(y, na.rm=TRUE);
  # Not needed anymore
  y <- NULL;
  verbose && exit(verbose);

  # Building fit
  fit <- list(mu=mu);

  verbose && exit(verbose);

  fit;
}, protected=TRUE)





setMethodS3("getNormalizeSignalsOne", "ScaleNormalization3", function(this, df, fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  df <- Arguments$getInstanceOf(df, "AffymetrixCelFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalizing one array according to model fit");
  verbose && cat(verbose, "Full name: ", getFullName(df));

  verbose && enter(verbose, "Getting algorithm parameters");
  params <- getParameters(this, expand=TRUE);
  cells <- params$cellsToUpdate;
  shift <- params$shift;
  targetAvg <- params$targetAvg;
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading signals");
#  y <- extractMatrix(df, units=cells, drop=TRUE, verbose=verbose);  # NOTE: 'units' :(
  y <- getData(df, indices=cells, fields="intensities", drop=TRUE, verbose=verbose);
  verbose && str(verbose, y);
  verbose && exit(verbose);

  # Shift?
  if (shift != 0) {
    verbose && enter(verbose, "Shifting signals");
    y <- y + shift;
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Rescaling");
  b <- targetAvg / fit$mu;
  verbose && printf(verbose, "Scale factor: %.2f\n", b);
  y <- b*y;
  verbose && cat(verbose, "Normalized signals:");
  verbose && str(verbose, y);
  # Sanity check
  yM <- median(y, na.rm=TRUE);
  verbose && printf(verbose, "Median after: %.2f\n", yM);
  verbose && exit(verbose);

  verbose && exit(verbose);

  y;
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
#   \item{...}{Not used.}
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
setMethodS3("process", "ScaleNormalization3", function(this, ..., skip=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Scale normalizing data set");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  dataSet <- getInputDataSet(this);

  # Get algorithm parameters
  verbose && enter(verbose, "Getting algorithm parameters");
  params <- getParameters(this, expand=TRUE);
  verbose && str(verbose, params);
  verbose && exit(verbose);

  # Get the output path
  outputPath <- getPath(this);

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing ", length(dataSet), " arrays");
  for (kk in seq_along(dataSet)) {
    verbose && enter(verbose, "Array #", kk);
    df <- dataSet[[kk]];
    verbose && print(verbose, df);

    filename <- basename(getPathname(df));
    filename <- gsub("[.]cel$", ".CEL", filename);  # Only output upper case!
    pathname <- Arguments$getWritablePathname(filename, path=outputPath);
    pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);

    # Already normalized?
    if (skip && isFile(pathname)) {
      verbose && cat(verbose, "Normalized data file already exists: ",
                                                                   pathname);
      verbose && exit(verbose);
      next;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit model
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Fitting model");
    fit <- fitOne(this, df=df, verbose=less(verbose));
    verbose && str(verbose, fit);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Normalize data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Normalizing for fitted effects");
    yN <- getNormalizeSignalsOne(this, df=df, fit=fit,
                                                  verbose=less(verbose));
    verbose && str(verbose, yN);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store normalized data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Writing normalized data");
    writeSignals(this, pathname=pathname, cells=params$cellsToUpdate,
                 intensities=yN, templateFile=df, verbose=less(verbose));
    # Not needed anymore
    yN <- NULL;
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    ## Create checksum file
    dfZ <- getChecksumFile(pathname)

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create result set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outputDataSet <- getOutputDataSet(this, force=TRUE, verbose=less(verbose));

  # Update the output data set
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);

  outputDataSet;
})


############################################################################
# HISTORY:
# 2008-07-25
# o Added fitOne() and getNormalizedSignalsOne().
# o Renamed getSubsetToAvg() to getSubsetToFit() to harmonize names.
# o Rewrote to make use of new ProbeLevelTransform2.R.
# o Extracted from old ScaleNormalization.R.
# 2008-02-19
# o Added support for "-X", "-Y" and "-XY" for argument 'subsetToAvg'.
# o Added support for constructor argument 'shift'.
# 2007-09-18
# o Updated all getCellIndices() to use 'unlist=TRUE'.
# 2007-04-16
# o Created.
############################################################################
