###########################################################################/**
# @RdocClass ReseqCrosstalkCalibration
#
# @title "The ReseqCrosstalkCalibration class"
#
# \description{
#  @classhierarchy
#
#  This class represents a calibration function that transforms the
#  probe-level signals such that the signals from the four nucleotides
#  (A, C, G, T) are orthogonal.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{An @see "AffymetrixCelSet".}
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{targetAvg}{The signal(s) that the average of the sum of the
#     probe quartets should have after calibration.}
#   \item{subsetToAvg}{The indices of the cells (taken as the intersect of
#     existing indices) used to calculate average in order to rescale to
#     the target average. If @NULL, all probes are considered.}
#   \item{mergeGroups}{A @logical ...}
#   \item{flavor}{A @character string specifying what algorithm is used
#     to fit the crosstalk calibration.}
#   \item{alpha, q, Q}{Additional arguments passed to
#    \code{fitMultiDimensionalCone()}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("ReseqCrosstalkCalibration", function(dataSet=NULL, ..., targetAvg=2200, subsetToAvg=NULL, mergeGroups=FALSE, flavor=c("sfit", "expectile"), alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extraTags <- NULL;

  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "AffymetrixCelSet");

    cdf <- getCdf(dataSet);

    # Argument 'targetAvg':
    if (!is.null(targetAvg)) {
      targetAvg <- Arguments$getDouble(targetAvg, range=c(0, Inf));
    }

    # Argument 'subsetToAvg':
    if (is.null(subsetToAvg)) {
    } else if (is.character(subsetToAvg)) {
      if (subsetToAvg %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'subsetToAvg': ", subsetToAvg);
      }
      extraTags <- c(extraTags, subsetToAvg=subsetToAvg);
    } else {
      subsetToAvg <- Arguments$getIndices(subsetToAvg, max=nbrOfCells(cdf));
      subsetToAvg <- unique(subsetToAvg);
      subsetToAvg <- sort(subsetToAvg);
    }
  }

  # Argument 'mergeGroups':
  if (!is.null(mergeGroups)) {
    mergeGroups <- Arguments$getLogical(mergeGroups);
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  extend(ProbeLevelTransform(dataSet=dataSet, ...), "ReseqCrosstalkCalibration",
    "cached.setsOfProbes" = NULL,
    "cached.subsetToAvgExpanded" = NULL,
    .targetAvg = targetAvg,
    .subsetToAvg = subsetToAvg,
    .mergeGroups = mergeGroups,
    .extraTags = extraTags,
    .flavor = flavor,
    .alpha = alpha,
    .q = q,
    .Q = Q
  )
})


setMethodS3("getAsteriskTags", "ReseqCrosstalkCalibration", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # 'mergeGroups' tag
  if (this$.mergeGroups) {
    mergeGroupsTag <- NULL;
  } else {
    mergeGroupsTag <- "byGroup";
  }
  tags <- c(tags, mergeGroupsTag);

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getSubsetToAvg", "ReseqCrosstalkCalibration", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  subsetToAvg <- this$.subsetToAvg;

  # Expand?
  if (is.character(subsetToAvg)) {
    if (subsetToAvg %in% c("-X", "-Y", "-XY")) {
      verbose && enter(verbose, "Identify subset of units from genome information");
      verbose && cat(verbose, "subsetToAvg: ", subsetToAvg);

      # Look up in cache
      subset <- this$.subsetToAvgExpanded;
      if (is.null(subset)) {
        dataSet <- getInputDataSet(this);
        cdf <- getCdf(dataSet);

        # Get the genome information (throws an exception if missing)
        gi <- getGenomeInformation(cdf);
        verbose && print(verbose, gi);

        # Identify units to be excluded
        if (subsetToAvg == "-X") {
          subset <- getUnitsOnChromosomes(gi, 23, .checkArgs=FALSE);
        } else if (subsetToAvg == "-Y") {
          subset <- getUnitsOnChromosomes(gi, 24, .checkArgs=FALSE);
        } else if (subsetToAvg == "-XY") {
          subset <- getUnitsOnChromosomes(gi, 23:24, .checkArgs=FALSE);
        }

        verbose && cat(verbose, "Units to exclude: ");
        verbose && str(verbose, subset);

        # Identify the cell indices for these units
        subset <- getCellIndices(cdf, units=subset,
                                 useNames=FALSE, unlist=TRUE);
        verbose && cat(verbose, "Cells to exclude: ");
        verbose && str(verbose, subset);

        # The cells to keep
        subset <- setdiff(1:nbrOfCells(cdf), subset);

        verbose && cat(verbose, "Cells to include: ");
        verbose && str(verbose, subset);

        # Store
        this$.subsetToAvgExpanded <- subset;
      }

      subsetToAvg <- subset;
      # Not needed anymore
      subset <- NULL;

      verbose && exit(verbose);
    }
  }

  subsetToAvg;
}, protected=TRUE);



setMethodS3("getParameters", "ReseqCrosstalkCalibration", function(this, expand=TRUE, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters", expand=expand);

  params <- c(params, list(
    targetAvg = this$.targetAvg,
    subsetToAvg = this$.subsetToAvg,
    mergeGroups = this$.mergeGroups,
    flavor = this$.flavor,
    alpha = this$.alpha,
    q = this$.q,
    Q = this$.Q
  ));

  # Expand?
  if (expand) {
    params$subsetToAvg <- getSubsetToAvg(this);
  }

  params;
}, protected=TRUE)



setMethodS3("rescaleByAll", "ReseqCrosstalkCalibration", function(this, yAll, params, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  targetAvg <- params$targetAvg;
  nt <- length(targetAvg);
  if (nt != 1) {
    throw("In order rescale towards a global target average ('rescaleBy' == \"all\"), argument 'targetAvg' must be a scalar: ", paste(targetAvg, collapse=","));
  }

  if (verbose) {
    enter(verbose, "Rescaling toward target average");
    cat(verbose, "Target average: ", targetAvg);
    if (!is.null(params$subsetToAvg)) {
      cat(verbose, "Using subset of cells for estimate of target average:");
      src <- attr(params$subsetToAvg, "src");
      if (is.null(src))
         src <- params$subsetToAvg;
      str(verbose, src);
    }
    cat(verbose, "yAll: ");
    str(verbose, yAll);
  }

  # Total number of values
  n0 <- length(yAll);

  # Average of *all* values
  yAvg0 <- median(yAll, na.rm=TRUE);

  if (!is.null(params$subsetToAvg)) {
    y <- yAll[params$subsetToAvg];
    n <- length(y);
    if (n == 0) {
      throw("Cannot rescale to target average. There are no cells to average over.");
    }
    yAvg <- median(y, na.rm=TRUE);
    verbose && printf(verbose, "yAvg (using %d/%.1f%% summed pairs): %.2f of %.2f (%.1f%%)\n", n, 100*n/n0, yAvg, yAvg0, 100*yAvg/yAvg0);
    # Not needed anymore
    y <- n <- NULL;
  } else {
    yAvg <- yAvg0;
    verbose && printf(verbose, "yAvg (100%%): %.2f\n", yAvg);
  }

  if (!is.finite(yAvg)) {
    throw("Cannot rescale to target average. Signal average is non-finite: ", yAvg);
  }

  b <- targetAvg/yAvg;
  verbose && printf(verbose, "Scale factor: %.2f\n", b);
  yAll <- b*yAll;
  fit <- list(all=list(b=b));
  verbose && exit(verbose);

  attr(yAll, "fit") <- fit;

  verbose && cat(verbose, "Rescaling parameter estimates:");
  verbose && str(verbose, fit);
  verbose && exit(verbose);

  yAll;
}, protected=TRUE)


setMethodS3("getSetsOfProbes", "ReseqCrosstalkCalibration", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  setsOfProbes <- this$.setsOfProbes;

  if (is.null(setsOfProbes)) {
    verbose && enter(verbose, "Identifying sets of cell tuples according to the CDF");
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    verbose && cat(verbose, "Chip type: ", getChipType(cdf));

    params <- getParameters(this);
    mergeGroups <- params$mergeGroups;
    verbose && cat(verbose, "Merging groups: ", mergeGroups);

    setsOfProbes <- getCellQuartets(cdf, mergeGroups=mergeGroups, verbose=verbose);

#    verbose && cat(verbose, "setsOfProbes:");
#    verbose && str(verbose, setsOfProbes);

    this$.setsOfProbes <- setsOfProbes;

    verbose && exit(verbose);
  }

  setsOfProbes;
}, protected=TRUE)



setMethodS3("fitOne", "ReseqCrosstalkCalibration", function(this, yAll, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting crosstalk model");

  params <- getParameters(this);

  verboseL <- (verbose && isVisible(verbose, -50));
  verbose && enter(verbose, "Fitting model unit by unit and group by group");
  cellQuartets <- getSetsOfProbes(this, verbose=less(verbose, 20));

  nbrOfUnits <- length(cellQuartets);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  fits <- vector("list", nbrOfUnits);
  names(fits) <- names(cellQuartets);

  for (uu in seq_len(nbrOfUnits)) {
    keyUU <- names(cellQuartets)[uu];
    verbose && enter(verbose, sprintf("Unit #%d ('%s') of %d",
                                           uu, keyUU, nbrOfUnits));
    cellsUU <- cellQuartets[[uu]];
    nbrOfGroups <- length(cellsUU);
    verbose && cat(verbose, "Number of groups: ", nbrOfGroups);

    fitsUU <- vector("list", nbrOfGroups);
    names(fitsUU) <- names(cellsUU);

    for (gg in seq_len(nbrOfGroups)) {
      keyGG <- names(cellsUU)[gg];
      verbose && enter(verbose, sprintf("Group #%d ('%s') of %d",
                                gg, keyGG, nbrOfGroups));
      cellsGG <- cellsUU[[gg]];

      # Extracting data
      y <- yAll[cellsGG];
      dim(y) <- dim(cellsGG);

      # Fitting crosstalk model
      y <- t(y);
      fitGG <- fitMultiDimensionalCone(y, flavor=params$flavor,
                      alpha=params$alpha, q=params$q, Q=params$Q,
                                                   verbose=verboseL);

      fitsUU[[gg]] <- fitGG;
      verbose && exit(verbose);
    } # for (gg ...)

    fits[[uu]] <- fitsUU;
    verbose && exit(verbose);
  } # for (uu ...)
  # Not needed anymore
  fitsUU <- fitGG <- keyGG <- keyUU <- y <- cellsGG <- cellsUU <- NULL;

  verbose && exit(verbose);

##    callHooks(sprintf("%s.onFitOne", hookName), df=df, yAll=yAll, fits=fits, ...);
  fits;
}, protected=TRUE)  # fitOne()


setMethodS3("calibrateOne", "ReseqCrosstalkCalibration", function(this, yAll, fits, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Backtransforming (calibrating)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calibrating data");

  # Get algorithm parameters
  params <- getParameters(this);
  targetAvg <- params$targetAvg;

  # Identify cell tuples
  cellQuartets <- getSetsOfProbes(this, verbose=less(verbose, 20));

  verbose && enter(verbose, "Calibrating unit by unit and group by group");
  nbrOfUnits <- length(cellQuartets);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);

  rescaleFits <- vector("list", nbrOfUnits);
  names(rescaleFits) <- names(cellQuartets);

  for (uu in seq_len(nbrOfUnits)) {
    keyUU <- names(cellQuartets)[uu];
    verbose && enter(verbose, sprintf("Unit #%d ('%s') of %d",
                                           uu, keyUU, nbrOfUnits));
    cellsUU <- cellQuartets[[uu]];
    fitsUU <- fits[[uu]];
    nbrOfGroups <- length(cellsUU);
    verbose && cat(verbose, "Number of groups: ", nbrOfGroups);

    rescaleFitsUU <- vector("list", nbrOfGroups);
    names(rescaleFitsUU) <- names(fitsUU);
    for (gg in seq_len(nbrOfGroups)) {
      keyGG <- names(cellsUU)[gg];
      verbose && enter(verbose, sprintf("Group #%d ('%s') of %d",
                                           gg, keyGG, nbrOfGroups));
      cellsGG <- cellsUU[[gg]];
      fitGG <- fitsUU[[gg]];

      # Extracting data
      y <- yAll[cellsGG];
      dim(y) <- dim(cellsGG);

      # Calibrating
      y <- t(y);
      yC <- backtransformMultiDimensionalCone(y, fit=fitGG);
      # Not needed anymore
      y <- fitGG <- NULL;
      yC <- t(yC);

      # Rescaling toward target average?
      if (!is.null(targetAvg)) {
        verbose && enter(verbose, "Rescaling toward target average");
        yC <- rescaleByAll(this, yAll=yC, params=params, verbose=less(verbose));
        fit <- attr(yC, "fit");
        fit$params <- NULL;
  ##      fit$paramsShort <- paramsShort;

        rescaleFitsUU[[gg]] <- fit;
        # Not needed anymore
        fit <- NULL;
        verbose && exit(verbose);
      }

      # Updating
      yAll[cellsGG] <- yC;
      # Not needed anymore
      yC <- NULL;

      verbose && exit(verbose);
    } # for (gg ...)

    rescaleFits[[uu]] <- rescaleFitsUU;

    # Not needed anymore
    cellsUU <- fitsUU <- rescaleFitsUU <- NULL;
    verbose && exit(verbose);
  } # for (uu ...)

  verbose && exit(verbose);

  # Update fits
  attr(yAll, "rescaleFits") <- rescaleFits;

  verbose && exit(verbose);

  yAll;
}, protected=TRUE) # calibrateOne()



###########################################################################/**
# @RdocMethod process
#
# @title "Calibrates the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already calibrated is re-calibrated,
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
setMethodS3("process", "ReseqCrosstalkCalibration", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calibrating data set for cross talk");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already calibrated");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  mergeGroups <- params$mergeGroups;

  # To be retrieved when needed.
  cellQuartets <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For hooks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hookName <- "process.ReseqCrosstalkCalibration";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Precalculate some model fit parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Compressing model parameter to a short format");
  paramsShort <- params;
  paramsShort$subsetToAvg <- NULL;
#  paramsShort$subsetToAvgIntervals <- seqToIntervals(params$subsetToAvg);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(ds);
  nbrOfArrays <- length(ds);
  verbose && enter(verbose, "Calibrating ", nbrOfArrays, " arrays");
  verbose && cat(verbose, "Path: ", outputPath);
  for (kk in seq_len(nbrOfArrays)) {
    df <- ds[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                              kk, getName(df), nbrOfArrays));

    fullname <- getFullName(df);
    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Already calibrated?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Calibrated data file already exists: ", pathname);
    } else {
      if (is.null(cellQuartets)) {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Identify the cell-index matrix
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        verbose && enter(verbose, "Identifying cell-index matrices");
        cellQuartets <- getCellQuartets(cdf, mergeGroups=mergeGroups, verbose=verbose);
        verbose && cat(verbose, "cellQuartets:");
#        verbose && str(verbose, cellQuartets);
#        dim <- dim(cellQuartets);
        gc <- gc();
        verbose && print(verbose, gc);
        verbose && exit(verbose);
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading all probe intensities");
      yAll <- getData(df, fields="intensities", ...)$intensities;
      verbose && str(verbose, yAll);
      verbose && exit(verbose);

      modelFit <- list(
        paramsShort=paramsShort
      );

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Fitting crosstalk model
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      fits <- fitOne(this, yAll=yAll, cellQuartets=cellQuartets, verbose=verbose);

      # Store crosstalk model fit(s)
      modelFit$rccFits <- fits;

#return(list(modelFit=modelFit, yAll=yAll));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Backtransforming (calibrating)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      yAll <- calibrateOne(this, yAll=yAll, fits=fits, verbose=verbose);
      # Not needed anymore
      fits <- NULL;

      # Store rescaling parameters
      modelFit$rescaleFits <- attr(yAll, "rescaleFits");

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store model fit
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store fit and parameters (in case someone are interested in looking
      # at them later; no promises of backward compatibility though).
      filename <- sprintf("%s,fit.RData", fullname);
      fitPathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      saveObject(modelFit, file=fitPathname);
      verbose && str(verbose, modelFit, level=-50);
      # Not needed anymore
      modelFit <- NULL;


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing calibrated data");

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- (force && isFile(pathname));
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

      createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      # Write calibrated data to file
      verbose2 <- -as.integer(verbose)-2;
      .updateCel(pathnameT, intensities=yAll, verbose=verbose2);

      # Not needed anymore
      yAll <- verbose2 <- NULL;

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose);

      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
    }

    # Test if calibrated data file can be retrieved
    dfC <- newInstance(df, pathname);

    ## Create checksum file
    dfCZ <- getChecksumFile(dfC)

    # Not needed anymore
    df <- NULL; # Not needed anymore

    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);

  # Garbage collect
  # Not needed anymore
  ds <- cellQuartets <- NULL;

  gc <- gc();
  verbose && print(verbose, gc);

  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);

  invisible(outputDataSet);
}) # process()





############################################################################
# HISTORY:
# 2008-08-29
# o Created fitOne() and calibrateOne().
# o Added option 'mergeGroups=FALSE'.
# o Updated according to new getCellQuartets() of AffymetrixCdfFile.
# 2008-08-10
# o Created from AllelicCrosstalkCalibration.R.
############################################################################
