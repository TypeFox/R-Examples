setConstructorS3("SpatialRowColumnNormalization", function(..., spar=c(0.7,0.7), blockSizes=c(20,20), maxIter=5) {
  # Argument 'spar':
  spar <- Arguments$getDoubles(spar, range=c(0,Inf));

  # Argument 'h':
  blockSizes <- Arguments$getIntegers(blockSizes, range=c(1,Inf));


  extend(ProbeLevelTransform(...), "SpatialRowColumnNormalization",
    .spar = spar,
    .blockSizes = blockSizes,
    .maxIter = maxIter
  );
})


setMethodS3("getParameters", "SpatialRowColumnNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    spar = getSpar(this),
    blockSizes = getBlockSizes(this),
    maxIter = this$.maxIter
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)


setMethodS3("getSpar", "SpatialRowColumnNormalization", function(this, ...) {
  this$.spar;
})

setMethodS3("getBlockSizes", "SpatialRowColumnNormalization", function(this, ...) {
  this$.blockSizes;
})


setMethodS3("process", "SpatialRowColumnNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalizing data set spatially in blocks of rows and columns");

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
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);

  # Get algorithm parameters (including the target distribution above)
  params <- getParameters(this);

  # Get the output path
  outputPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(ds);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);
  dataFiles <- list();
  cells <- NULL;
  for (kk in seq_len(nbrOfArrays)) {
    df <- ds[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                              kk, getName(df), nbrOfArrays));

    fullname <- getFullName(df);
    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Already normalized?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Normalized data file already exists: ", pathname);
    } else {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading all probe intensities");
      y <- readRawDataRectangle(df, field="intensities", drop=TRUE, ...);
      dim <- dim(y);
      verbose && str(verbose, y);
      verbose && exit(verbose);

      verbose && enter(verbose, "Transforming to the log-scale");
      y <- log2(y);
      verbose && str(verbose, y);
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Normalizing log-ratios data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Normalizing rows and columns in blocks");
      fit <- fitSplineBlockPolish(y, blockSizes=params$blockSizes,
                              spar=params$spar, maxIter=params$maxIter, ...);
      verbose && str(verbose, fit);
      y <- residuals(fit);
      # Not needed anymore
      fit <- NULL;
      verbose && exit(verbose);

      verbose && enter(verbose, "Back-transforming to intensity scale");
      y <- y + 12;
      y <- 2^y;
      y <- as.vector(y);
      verbose && str(verbose, y);
      verbose && exit(verbose);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing normalized data");

      if (is.null(cells)) {
        cells <- matrix(seq_along(y), nrow=dim[1], ncol=dim[2], byrow=TRUE);
        verbose && str(verbose, cells);
      }

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- (force && isFile(pathname));
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

      createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      # Write calibrated data to file
      verbose2 <- -as.integer(verbose)-2;
      .updateCel(pathnameT, indices=cells, intensities=y, verbose=verbose2);
      # Not needed anymore
      y <- verbose2 <- NULL;

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose);

      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);
    } # if-else

    # Retrieving normalized data file
    dfN <- newInstance(df, pathname);

    # CDF inheritance
    setCdf(dfN, cdf);

    # Record
    dataFiles[[kk]] <- dfN;

    # Not needed anymore
    df <- dfN <- NULL;
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Garbage collect
  # Not needed anymore
  dataFiles <- ds <- cells <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  # Update the output data set
  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);

  invisible(outputDataSet);
})


############################################################################
# HISTORY:
# 2008-04-02
# o Note: Normalizing log-ratios and then transforming back to chip effects
#   might not do.  See my PPT notes from today.
# 2008-03-19
# o Created.
############################################################################
