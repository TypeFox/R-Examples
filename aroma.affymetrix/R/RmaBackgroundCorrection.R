###########################################################################/**
# @RdocClass RmaBackgroundCorrection
#
# @title "The RmaBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents the RMA background adjustment function.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "BackgroundCorrection".}
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
#   The fitting algorithm of the RMA background correction model may not
#   converge if there too many small and discrete signals.  To overcome
#   this problem, a small amount of noise may be added to the signals
#   before fitting the model.  This is an ad hoc solution that seems to
#   work.
#   However, add Gaussian noise may generate non-positive signals.
# }
#
# \details{
#   Internally @see "affy::bg.adjust" is used to background correct the
#   probe signals.  The default is to background correct PM signals only.
# }
#
# @author "KS, HB"
#*/###########################################################################
setConstructorS3("RmaBackgroundCorrection", function(..., addJitter=FALSE, jitterSd=0.2, seed=6022007) {
  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  extend(BackgroundCorrection(..., typesToUpdate="pm"),
    "RmaBackgroundCorrection",
    .addJitter=addJitter,
    .jitterSd=jitterSd,
    .seed=seed
  );
})


setMethodS3("getParameters", "RmaBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  pmOnly <- (this$.typesToUpdate == "pm");

  # Get parameters of this class
  params2 <- list(
    addJitter = this$.addJitter,
    jitterSd = this$.jitterSd,
    seed = this$.seed,
    pmonly = pmOnly
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)



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
#  Returns the output data set.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "RmaBackgroundCorrection", function(this, ..., force=FALSE, verbose=FALSE) {
  # Load required packages
  requireNamespace("affy") || throw("Package not loaded: affy")
  bg.adjust <- affy::bg.adjust

  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Background correcting data set")

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already background corrected")
    verbose && exit(verbose)
    outputDataSet <- getOutputDataSet(this)
    return(outputDataSet)
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this)

  # Get the output path
  outputPath <- getPath(this)

  cdf <- getCdf(ds)

  # Get algorithm parameters (including the target distribution)
  params <- getParameters(this)
  # 'subsetToUpdate' is not used and 'typesToUpdate' are used via 'pmonly'
  pmonly <- params$pmonly
  pmCells <- NULL
  pmJitter <- NULL
  addJitter <- params$addJitter
  jitterSd <- params$jitterSd
  # Not needed anymore
  params <- NULL


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Background correct
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(ds)
  verbose && cat(verbose, "Number of arrays: ", nbrOfArrays)

  res <- listenv()

  for (ii in seq_along(ds)) {
    df <- ds[[ii]]
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", ii, getName(df), nbrOfArrays))

    filename <- basename(getPathname(df))
    filename <- gsub("[.]cel$", ".CEL", filename)
    pathname <- Arguments$getWritablePathname(filename, path=outputPath,
                                                        mustNotExist=FALSE)
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already processed. Skipping.")

      ## Assert validity
      dfD <- newInstance(df, pathname)
      setCdf(dfD, cdf)

      res[[ii]] <- pathname
      verbose && exit(verbose)
      next
    }

    if (pmonly && is.null(pmCells)) {
      verbose && enter(verbose, "Identifying PM-only probes (pmonly=TRUE)")
      verbose && print(verbose, cdf)
      chipType <- getChipType(cdf)
      key <- list(method="bgAdjustRma", class=class(cdf)[1], chipType=chipType)
      dirs <- c("aroma.affymetrix", chipType)
      pmCells <- loadCache(key=key, dirs=dirs)
      if (is.null(pmCells)) {
        indices <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE)
        pmCells <- indices[isPm(cdf)]
        saveCache(pmCells, key=key, dirs=dirs)
      }
      nbrOfPMs <- length(pmCells)
      verbose && exit(verbose)
    } else {
      nbrOfPMs <- nbrOfCells(cdf)
    }

    if (addJitter && is.null(pmJitter)) {
      ## Use a temporary random seed?
      seed <- params$seed
      if (!is.null(seed)) {
        randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
        on.exit(randomSeed("reset"), add=TRUE)
        verbose && printf(verbose, "Random seed temporarily set (seed=%d, kind=\"L'Ecuyer-CMRG\")\n", seed)
      }
      pmJitter <- rnorm(nbrOfPMs, mean=0, sd=jitterSd)
    }

    res[[ii]] %<=% {
      verbose && enter(verbose, "Obtaining signals")
      pm <- readRawData(df, indices=pmCells, "intensities")$intensities
      if (addJitter) pm <- pm + pmJitter
      clearCache(df)
      verbose && exit(verbose)

      # adjust background - use original affy functions to avoid errors from
      # re-coding
      verbose && enter(verbose, "Applying normal+exponential signal model")
      pm <- bg.adjust(pm)  # From package 'affy' (without a namespace)
      verbose && exit(verbose)

      # update the PM

      # Write adjusted data to file
      verbose && enter(verbose, "Writing adjusted probe signals")

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- (force && isFile(pathname))
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose)

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing")
      createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose))
      verbose && exit(verbose)

      verbose && enter(verbose, "Writing adjusted intensities")
      .updateCel(pathnameT, indices=pmCells, intensities=pm)
      verbose && exit(verbose)

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose)

      ## Create checksum file
      dfZ <- getChecksumFile(pathname)

      verbose && exit(verbose)

      # Not needed anymore
      pm <- NULL

      # Garbage collection
      gc <- gc()
      verbose && print(verbose, gc)

      pathname
    } ## %<=%

    verbose && exit(verbose)
  } # for (ii ...)

  ## Not needed anymore
  pmCells <- NULL

  ## Resolve futures
  res <- as.list(res)
  res <- NULL

  ## Garbage collect
  gc <- gc()
  verbose && print(verbose, gc)

  # Get the output data set
  outputDataSet <- getOutputDataSet(this)

  verbose && exit(verbose)

  outputDataSet
})


############################################################################
# HISTORY:
# 2015-11-19
# o CLEANUP: Migrated bgAdjustRma() for AffymetrixCelFile into process().
# o CLEANUP: Using readRawData() instead of getData().
# 2012-11-20
# o CLEANUP: process() for RmaBackgroundCorrection now processes
#   each file by itself, i.e. it no longer calls bgAdjustRma() for
#   AffymetrixCelSet (which has been removed).
# 2007-06-30
# o Added Rdoc comments about jitter.
# 2007-05-26
# o Updated the Rdocs.
# 2007-03-21
# o Created.
############################################################################
