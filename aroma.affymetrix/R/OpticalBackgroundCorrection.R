###########################################################################/**
# @RdocClass OpticalBackgroundCorrection
#
# @title "The OpticalBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents "optical" background adjustment.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{minimum}{The minimum signal allowed after adjustment.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS"
#*/###########################################################################
setConstructorS3("OpticalBackgroundCorrection", function(..., minimum=1) {
  extend(BackgroundCorrection(..., typesToUpdate="pmmm"),
    "OpticalBackgroundCorrection",
    .minimum = minimum
  )
})


setMethodS3("getParameters", "OpticalBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    minimum = this$.minimum
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
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "OpticalBackgroundCorrection", function(this, ..., force=FALSE, verbose=FALSE) {

  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Background correcting data set")

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already background corrected for \"optical\" effects")
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
  subsetToUpdate <- params$subsetToUpdate
  typesToUpdate <- params$typesToUpdate
  minimum <- params$minimum
  # Not needed anymore
  params <- NULL


  hasSubsetToUpdate <- getFraction <- FALSE
  nbrOfCells <- nbrOfCells(cdf)

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
    filename <- gsub("[.]cel$", ".CEL", filename);  # Only output upper case!
    pathname <- Arguments$getWritablePathname(filename, path=outputPath,
                                              mustNotExist=FALSE)
    pathname <- AffymetrixFile$renameToUpperCaseExt(pathname)

    # Already processed?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already processed. Skipping.")
      # Assert valid file
      dfOut <- newInstance(df, pathname)
      setCdf(dfOut, cdf)
      res[[ii]] <- pathname
      verbose && exit(verbose)
      next
    }


    # Identifying the cells to be updated?
    if (!hasSubsetToUpdate) {
      verbose && enter(verbose, "Identifying cells to be updated")
      subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate,
                                                          types=typesToUpdate)
      verbose && cat(verbose, "Number of cells: ", length(subsetToUpdate))
      verbose && str(verbose, subsetToUpdate)
      hasSubsetToUpdate <- TRUE
      verbose && exit(verbose)
    }


    res[[ii]] %<=% {
      # Get all probe signals
      verbose && enter(verbose, "Reading probe intensities")
      x <- readRawData(df, fields="intensities", verbose=less(verbose,2))
      x <- x$intensities
      verbose && exit(verbose)

      # Subtract optical background from selected probes
      verbose && enter(verbose, "Adjusting background for optical effect")
      arrayMinimum <- min(x[subsetToUpdate], na.rm=TRUE)
      verbose && printf(verbose, "Array minimum: %.2f\n", arrayMinimum)
      xdiff <- (arrayMinimum - minimum)
      verbose && printf(verbose, "Correction: -(%.2f-%.2f) = %+.2f\n",
                                             arrayMinimum, minimum, -xdiff)
      x[subsetToUpdate] <- x[subsetToUpdate] - xdiff
      verbose && exit(verbose)

      # Write adjusted data to file
      verbose && enter(verbose, "Writing adjusted probe signals")

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- (force && isFile(pathname));
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose)

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing")
      createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose))
      verbose && exit(verbose)

      verbose && enter(verbose, "Writing adjusted intensities")
      .updateCel(pathnameT, intensities=x)
      verbose && exit(verbose)

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose)

      dfOut <- newInstance(df, pathname)
      setCdf(dfOut, cdf)

      ## Create checksum file
      dfZ <- getChecksumFile(dfOut)

      verbose && exit(verbose)

      pathname
    } ## %<=%

    verbose && exit(verbose)
  } # for (ii ...)

  ## Not needed anymore
  subsetToUpdate <- NULL

  ## Resolve futures
  res <- as.list(res)
  res <- NULL

  ## Garbage collect
  gc <- gc()
  verbose && print(verbose, gc)

  ## Get the output data set
  outputDataSet <- getOutputDataSet(this)

  verbose && exit(verbose)

  outputDataSet
})



############################################################################
# HISTORY:
# 2012-11-20
# o CLEANUP: process() for OpticalBackgroundCorrection now processes
#   each file by itself, i.e. it no longer calls bgAdjustOptical() for
#   AffymetrixCelSet (which has been removed).
# 2007-08-24
# o BUG FIX: Forgot to pass argument '.deprecated=FALSE' to bgAdjustGcrma()
#   because the latter is deprecated at the user-API level.
# 2007-03-22
# o Created.
############################################################################
