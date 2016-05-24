###########################################################################/**
# @RdocClass LimmaBackgroundCorrection
#
# @title "The LimmaBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents the various "background" correction methods
#  implemented in the \pkg{limma} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "BackgroundCorrection".}
#   \item{args}{A @list of additional arguments passed to the
#     correction algorithm.}
#   \item{addJitter}{If @TRUE, Zero-mean gaussian noise is added to the
#     signals before being background corrected.}
#   \item{jitterSd}{Standard deviation of the jitter noise added.}
#   \item{seed}{An (optional) @integer specifying a temporary random seed
#     to be used for generating the (optional) jitter.  The random seed
#     is set to its original state when done.  If @NULL, it is not set.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Jitter noise}{
#   The fitting algorithm of the normal+exponentital background correction
#   model may not converge if there too many small and discrete signals.
#   To overcome this problem, a small amount of noise may be added to the
#   signals before fitting the model.  This is an ad hoc solution that
#   seems to work.
#   However, adding Gaussian noise may generate non-positive signals.
# }
#
# \details{
#   By default, only PM signals are background corrected and MMs are
#   left unchanged.
# }
#
# \author{Henrik Bengtsson.
#         Adopted from RmaBackgroundCorrection by Ken Simpson.}
#
# \seealso{
#   Internally, @see "limma::backgroundCorrect" is used.
# }
#
#*/###########################################################################
setConstructorS3("LimmaBackgroundCorrection", function(..., args=NULL, addJitter=FALSE, jitterSd=0.2, seed=6022007) {
  # Argument 'args':
  if (!is.null(args)) {
    if (!is.list(args)) {
      throw("Argument 'args' is not a list: ", class(args)[1]);
    }
  }

  # Argument 'addJitter':
  addJitter <- Arguments$getLogical(addJitter);

  # Argument 'jitterSd':
  jitterSd <- Arguments$getDouble(jitterSd);

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  extend(BackgroundCorrection(..., typesToUpdate="pm"), "LimmaBackgroundCorrection",
    .args = args,
    .addJitter = addJitter,
    .jitterSd = jitterSd,
    .seed = seed
  );
})

setMethodS3("getAsteriskTags", "LimmaBackgroundCorrection", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Extra tags?
  params <- getParameters(this);
  args <- params$args;
  tags <- c(tags, args$method, args$normexp.method);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getParameters", "LimmaBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  pmOnly <- (this$.typesToUpdate == "pm");

  # Get parameters of this class
  params2 <- list(
    addJitter = this$.addJitter,
    jitterSd = this$.jitterSd,
    pmOnly = pmOnly,
    seed = this$.seed
  );

  # Algorithm parameters
  params$args <- this$.args;

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)


setMethodS3("getSubsetToUpdate0", "LimmaBackgroundCorrection", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  cells <- this$.cellsToUpdate0;

  if (is.null(cells)) {
    # Get algorithm parameters
    params <- getParameters(this);

    if (params$pmOnly) {
      verbose && enter(verbose, "Retrieving PM-only CDF indices");
      ds <- getInputDataSet(this);
      cdf <- getCdf(ds);
      chipType <- getChipType(cdf);
      key <- list(method="getPmCellIndices", class=class(cdf)[1], chipType=chipType);
      dirs <- c("aroma.affymetrix", chipType);
      cells <- loadCache(key=key, dirs=dirs);
      if (is.null(cells)) {
        indices <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE);
        cells <- indices[isPm(cdf)];
        # Not needed anymore
        indices <- NULL;
        saveCache(cells, key=key, dirs=dirs);
      }
      verbose && exit(verbose);
    }

    this$.cellsToUpdate0 <- cells;
  }

  cells;
}, private=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Performs background correction"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
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
setMethodS3("process", "LimmaBackgroundCorrection", function(this, ..., force=FALSE, verbose=FALSE) {

  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Background correcting data set");

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already background corrected");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get the CDF
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get the output path
  outputPath <- getPath(this);

  # The cells to be updated (defaults to all, but will be overloaded
  # below in the typical settings)
  cells <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Try to load the required packaged
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  requireNamespace("limma") || throw("Package not loaded: limma")
  backgroundCorrect <- limma::backgroundCorrect

  ## In case random jitter will be used
  jitter <- NULL


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # apply normal+exponential model to each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(ds);
  verbose && enter(verbose, "Adjusting ", nbrOfArrays, " arrays");

  res <- listenv()

  for (kk in seq_along(ds)) {
    verbose && enter(verbose, sprintf("Array #%d of %d", kk, nbrOfArrays));
    df <- ds[[kk]];
    verbose && print(verbose, df);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Generating output pathname
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fullname <- getFullName(df);
    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Nothing to do? Already corrected?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Output data file already exists: ", pathname);
      verbose && exit(verbose);
      res[[kk]] <- pathname
      next;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify the indices for cells to be corrected
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(cells)) {
      cells <- getSubsetToUpdate0(this, verbose=less(verbose, 10));
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Generate random jitter?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (params$addJitter && is.null(jitter)) {
      ## Use a temporary random seed?
      seed <- params$seed
      if (!is.null(seed)) {
        randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
        on.exit(randomSeed("reset"), add=TRUE)
        verbose && printf(verbose, "Random seed temporarily set (seed=%d, kind=\"L'Ecuyer-CMRG\")\n", seed)
      }
      jitter <- rnorm(length(cells), mean=0, sd=params$jitterSd);
    }


    res[[kk]] %<=% {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Extracting data");
      y <- extractMatrix(df, cells=cells, drop=TRUE);
      verbose && str(verbose, y);
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Add jitter?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (params$addJitter) {
        y <- y + jitter;
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Correct data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Calling backgroundCorrect() of limma");
      args <- c(list(y), params$args);
      verbose && str(verbose, args);
      y <- do.call(backgroundCorrect, args=args);
      verbose && str(verbose, y);
      y <- y[,1,drop=TRUE];
      verbose && str(verbose, y);
      verbose && exit(verbose);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing corrected data");

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- (force && isFile(pathname));
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      # Write calibrated data to file
      verbose2 <- -as.integer(verbose)-2;
      .updateCel(pathnameT, indices=cells, intensities=y, verbose=verbose2);
      # Not needed anymore
      y <- verbose2 <- NULL;

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose);

      ## Create checksum file
      dfZ <- getChecksumFile(pathname)

      verbose && exit(verbose);

      pathname
    } ## %<=%

    # Not needed anymore
    df <- NULL;
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  ## Not needed anymore
  cells <- jitter <- NULL

  ## Resolve futures
  res <- as.list(res)
  res <- NULL

  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);

  invisible(outputDataSet);
})


############################################################################
# HISTORY:
# 2009-04-16
# o Made this a limma-only class. Removed the 'flavor' argument.
# 2009-04-09
# o Added redundancy test for LimmaBackgroundCorrection.
# 2009-04-06
# o Verified that new LimmaBackgroundCorrection and old
#   RmaBackgroundCorrection gives identical results.
# o Added support for 'affy' and 'limma' flavor.
# o Added support to specify and pass any algorithm parameters.
# o Created from RmaBackgroundCorrection.R.
############################################################################
