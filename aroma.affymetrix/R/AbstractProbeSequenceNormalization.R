###########################################################################/**
# @RdocClass AbstractProbeSequenceNormalization
#
# @title "The AbstractProbeSequenceNormalization class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a normalization method that corrects for
#  systematic effects in the probe intensities due to differences in
#  probe sequences.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform3".}
#   \item{target}{A @character string specifying type of "target" used.
#     If \code{"zero"}, all arrays are normalized to have no effects.
#     If @NULL, all arrays a normalized to have the same effect as
#     the average array has.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires that an @see "aroma.core::AromaCellSequenceFile" is
#   available for the chip type.
# }
#
# @author
#*/###########################################################################
setConstructorS3("AbstractProbeSequenceNormalization", function(..., target=NULL) {
  # Argument 'target':
  if (is.null(target)) {
  } else {
    target <- Arguments$getCharacter(target, length=c(1,1));
    knowTargets <- c("zero");
    if (!target %in% knowTargets) {
      throw("The value of argument 'target' is unknown: ", target);
    }
  }

  extend(ProbeLevelTransform3(...), "AbstractProbeSequenceNormalization",
    .target = target
  );
}, abstract=TRUE)



setMethodS3("getParameters", "AbstractProbeSequenceNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  params <- c(params, list(
    target = this$.target
  ));

  params;
}, protected=TRUE)



setMethodS3("getTargetFile", "AbstractProbeSequenceNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  dataSet <- getInputDataSet(this);
  dfR <- getAverageFile(dataSet, verbose=less(verbose, 25));

  dfR;
})


setMethodS3("getAromaCellSequenceFile", "AbstractProbeSequenceNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting AromaCellSequenceFile");

  aps <- this$.aps;

  if (force || is.null(aps)) {
    verbose && enter(verbose, "Locating");

    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf, fullname=FALSE);
    nbrOfCells <- nbrOfCells(cdf);
    # Not needed anymore
    dataSet <- cdf <- NULL;

    verbose && cat(verbose, "Chip type: ", chipType);
    verbose && cat(verbose, "Number of cells: ", nbrOfCells);

    aps <- AromaCellSequenceFile$byChipType(chipType,
                            nbrOfCells=nbrOfCells, ..., verbose=verbose);

    verbose && exit(verbose);

    this$.aps <- aps;
  }

  verbose && print(verbose, aps);
  verbose && exit(verbose);

  aps;
}, protected=TRUE)



setMethodS3("indexOfMissingSequences", "AbstractProbeSequenceNormalization", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Identifying probes with missing sequences");
  # Locate AromaCellSequenceFile holding probe sequences
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));

  idxs <- isMissing(acs, verbose=less(verbose, 5));
  idxs <- which(idxs);
  verbose && cat(verbose, "Cells with unknown sequences:");
  verbose && str(verbose, idxs);
  # Not needed anymore
  acs <- NULL;
  verbose && exit(verbose);

  idxs;
}, protected=TRUE);


setMethodS3("fitOne", "AbstractProbeSequenceNormalization", abstract=TRUE, protected=TRUE);

setMethodS3("predictOne", "AbstractProbeSequenceNormalization", abstract=TRUE, protected=TRUE);


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
#   \item{ram}{A positive @double scale factor specifying how much more
#     memory to use relative to the default.}
#   \item{force}{If @TRUE, data already normalized is re-normalized,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "AbstractProbeSequenceNormalization", function(this, ..., ram=NULL, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readSeqs <- function(this, cells, ...) {
    verbose && enter(verbose, "Reading probe sequences");
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);

    acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
    seqs <- readSequenceMatrix(acs, cells=cells, what="raw",
                                                  verbose=less(verbose, 5));
    # Not needed anymore
    acs <- NULL;
    verbose && exit(verbose);
    seqs;
  } # readSeqs()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Normalizing data set for probe-sequence effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    ## FIXME: outputDataSetZ <- getChecksumFileSet(outputDataSet)
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this, expand=TRUE, verbose=less(verbose, 5));

  # Get (and create) the output path
  outputPath <- getPath(this);

  # Get target
  target <- params$target;

  # Get shift
  shift <- params$shift;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(ds);
  df <- getOneFile(ds);
  nbrOfCells <- nbrOfCells(df);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);
  # Not needed anymore
  df <- NULL;

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  paramsShort <- NULL;
  muT <- NULL;
  seqs <- NULL;

  res <- listenv()

  for (kk in seq_len(nbrOfArrays)) {
    df <- ds[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                              kk, getName(df), nbrOfArrays));

    fullname <- getFullName(df);
    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Already calibrated?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Normalized data file already exists: ", pathname)
      res[[kk]] <- pathname
      verbose && exit(verbose)
      next
    }


    if (is.null(paramsShort)) {
      # Precalculate some model fit parameters
      verbose && enter(verbose, "Compressing model parameter to a short format");
      paramsShort <- params;
      paramsShort$cellsToFit <- NULL;
      paramsShort$cellsToUpdate <- NULL;
##      paramsShort$unitsToFitIntervals <- seqToIntervals(paramsShort$unitsToFit);
##      paramsShort$unitsToUpdateIntervals <- seqToIntervals(paramsShort$unitsToUpdate);
      verbose && exit(verbose);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Setting up model fit parameters
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    modelFit <- list(
      paramsShort = paramsShort
    )



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Phase 0: Fit probe-sequence effect for target?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(target) && is.null(muT)) {
      verbose && enter(verbose, "Modelling effects of target array");

      modelFitT <- list(
        paramsShort = paramsShort
      );

      dfT <- getTargetFile(this, verbose=less(verbose, 5));
      fullnameT <- getFullName(dfT);
      filename <- sprintf("%s,fit.RData", fullnameT);
      fitPathname <- Arguments$getWritablePathname(filename,
                                                     path=outputPath, ...);
      if (isFile(fitPathname)) {
        verbose && enter(verbose, "Loading already fitted probe-sequence effects for target");
        verbose && cat(verbose, "Pathname: ", fitPathname);
        modelFitT <- loadObject(file=fitPathname);
        fitT <- modelFitT$fit;
        verbose && exit(verbose);
      } else {
        verbose && enter(verbose, "Estimating probe-sequence effects for target");
        fitT <- fitOne(this, df=dfT, params=params, ram=ram, verbose=less(verbose, 5));
        verbose && print(verbose, fitT);
        modelFitT$fit <- fitT;

        verbose && enter(verbose, "Saving model fit");
        # Store fit and parameters (in case someone are interested in looking
        # at them later; no promises of backward compatibility though).
        saveObject(modelFitT, file=fitPathname);
        verbose && exit(verbose);

        verbose && exit(verbose);
      }

      verbose && str(verbose, modelFitT, level=-50);
      # Not needed anymore
      dfT <- modelFitT <- NULL;
      verbose && exit(verbose);


      verbose && enter(verbose, "Predicting probe affinities");
      if (is.null(seqs)) {
        seqs <- readSeqs(this, cells=params$cellsToUpdate);
      }
      muT <- predictOne(this, fit=fitT, params=params, seqs=seqs, verbose=less(verbose, 5));
      # Not needed anymore
      fitT <- NULL;
      verbose && cat(verbose, "muT:");
      verbose && str(verbose, muT);
      verbose && summary(verbose, muT);
      verbose && exit(verbose);

      # Garbage collection
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
    } ## if (is.null(target) && is.null(muT))

    if (is.null(seqs)) {
      seqs <- readSeqs(this, cells=params$cellsToUpdate);
    }


    res[[kk]] %<=% {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Future related
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Help identifying some globals (required)
      modelFit

      ## Prevent params from being exported as a global (optional)
      params <- getParameters(this, expand=TRUE)


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Phase I: Fit probe-sequence effect for the current array
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Fitting model for current array");
      fit <- fitOne(this, df=df, params=params, ram=ram, verbose=less(verbose, 5));
      verbose && print(verbose, fit);
      modelFit$fit <- fit;
      verbose && exit(verbose);

      verbose && enter(verbose, "Saving model fit");
      # Store fit and parameters (in case someone are interested in looking
      # at them later; no promises of backward compatibility though).
      filename <- sprintf("%s,fit.RData", fullname);
      fitPathname <- Arguments$getWritablePathname(filename,
                                                      path=outputPath, ...);
      saveObject(modelFit, file=fitPathname);
      verbose && str(verbose, modelFit, level=-50);
      # Not needed anymore
      modelFit <- NULL;
      verbose && exit(verbose);

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);



      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Phase II: Normalize current array toward target
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading probe signals");
      y <- extractMatrix(df, cells=params$cellsToUpdate, drop=TRUE);

      # Shift signals?
      if (shift != 0) {
        y <- y + shift;
        verbose && cat(verbose, "Shifted probe signals: ", shift);
      }

      verbose && str(verbose, y);
      verbose && summary(verbose, y);
      verbose && exit(verbose);

      verbose && enter(verbose, "Predicting mean log2 probe signals");
      mu <- predictOne(this, fit=fit, params=params, seqs=seqs, verbose=less(verbose, 5));
      # Not needed anymore
      fit <- seqs <- NULL
      verbose && cat(verbose, "mu:");
      verbose && str(verbose, mu);
      verbose && summary(verbose, mu);

      verbose && exit(verbose);


      verbose && enter(verbose, "Discrepancy scale factors towards target");
      verbose && cat(verbose, "Target: ", target);
      if (is.null(target)) {
        rho <- (muT-mu);
      } else if (target == "zero") {
        rho <- -mu;
      }
      # Not needed anymore
      mu <- NULL;
      verbose && summary(verbose, rho);
      rho <- 2^rho;
      verbose && summary(verbose, rho);

      # Update only subset with "finite" corrections
      keep <- which(is.finite(rho));
      rho <- rho[keep];
      y <- y[keep];
      cellsToUpdateKK <- params$cellsToUpdate[keep];
      # Not needed anymore
      keep <- NULL;
      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      verbose && enter(verbose, "Normalizing probe signals");
      y <- rho * y;
      # Not needed anymore
      rho <- NULL;
      verbose && str(verbose, y);
      verbose && summary(verbose, y);
      verbose && exit(verbose);

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing normalized data");

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- (force && isFile(pathname));
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      # Write calibrated data to file
      verbose2 <- -as.integer(verbose)-2;
      .updateCel(pathnameT, indices=cellsToUpdateKK, intensities=y, verbose=verbose2);
      # Not needed anymore
      y <- cellsToUpdateKK <- verbose2 <- NULL;
      gc <- gc();
      verbose && print(verbose, gc);

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose);

      verbose && exit(verbose);

      # Validating by retrieving calibrated data file
      dfC <- newInstance(df, pathname);

      ## Generate checksum file
      dfCZ <- getChecksumFile(dfC)

      # Not needed anymore
      dfC <- dfCZ <- NULL

      # Garbage collection
      gc <- gc();
      verbose && print(verbose, gc);

      pathname
    } ## %<=%

    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);

  ## Not needed anymore
  ds <- seqs <- muT <- NULL;

  ## Resolve futures
  res <- as.list(res)
  res <- NULL

  ## Garbage collect
#  clearCache(this);
  gc <- gc();
  verbose && print(verbose, gc);

  outputDataSet <- getOutputDataSet(this, force=TRUE)
  ## FIXME: outputDataSetZ <- getChecksumFileSet(outputDataSet)

  verbose && exit(verbose);

  invisible(outputDataSet);
})



############################################################################
# HISTORY:
# 2011-02-15
# o Added forgot space is two verbose statements.
# 2010-02-16
# o Added verbose output to getAromaCellSequenceFile().
## 2010-02-15
## o MEMORY OPTIMIZATION: Now process() of AbstractProbeSequenceNormalization
##   clears the in-memory cache when finished.
# 2009-07-08
# o ROBUSTNESS: Updated process() of AbstractProbeSequenceNormalization to
#   write to a tempory file which is the renamed.  This will lower the risk
#   for generating corrupt files in case of interrupts.
# 2008-12-03
# o SPEED UP: Now the "expanded" algorithm parameters ('params') are passed
#   to fitOne() and predictOne().  It is up to the implementation of these
#   two to either use them or not.
# 2008-11-29
# o Added argument 'ram' to process().
# 2008-08-05
# o Added support for specifying the type of target effects for any
#   AbstractProbeSequenceNormalization method.
# 2008-07-28
# o Now the model fit for the target is also saved to file.
# 2008-07-21
# o The process() is rather generic. Subclasses have to implement fitOne()
#   and predictOne().
# o Created from BaseCountNormalization.R.
############################################################################
