###########################################################################/**
# @RdocClass UnitTypeScaleNormalization
#
# @title "The UnitTypeScaleNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization function that transforms the
#  probe signals such that each unit type gets the same average.
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
setConstructorS3("UnitTypeScaleNormalization", function(..., targetAvg=4400) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  targetAvg <- Arguments$getDouble(targetAvg, range=c(1,Inf));

  extend(ProbeLevelTransform3(...), "UnitTypeScaleNormalization",
    .targetAvg = targetAvg
  )
})



setMethodS3("getParameters", "UnitTypeScaleNormalization", function(this, ...) {
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



setMethodS3("fitOne", "UnitTypeScaleNormalization", function(this, df, ..., verbose=FALSE) {
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
  units <- params$unitsToFit;
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  shift <- params$shift;
  stratifyBy <- params$typesToFit;
  verbose && exit(verbose);

  verbose && enter(verbose, "Getting unit types");
  cdf <- getCdf(df);
  unitTypes <- getUnitTypes(cdf, units=units, verbose=verbose);
  uniqueUnitTypes <- sort(unique(unitTypes));
  knownUnitTypes <- attr(unitTypes, "types");
  verbose && exit(verbose);

  fit <- list();
  for (kk in seq_along(uniqueUnitTypes)) {
    unitType <- uniqueUnitTypes[kk];
    key <- names(knownUnitTypes)[unitType == knownUnitTypes];
    names(unitType) <- key;
    verbose && enter(verbose, sprintf("Unit type %d (%d; '%s') of %d",
                    kk, unitType, key, length(uniqueUnitTypes)));
    unitsKK <- units[unitTypes == unitType];
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, unitsKK);
    nbrOfUnits <- length(unitsKK);

    verbose && enter(verbose, "Identifying cell indices");
    cells <- getCellsToInternal(this, units=unitsKK, stratifyBy=stratifyBy, verbose=verbose);
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);
    verbose && exit(verbose);
    nbrOfCells <- length(cells);

    verbose && enter(verbose, "Reading signals");
    y <- getData(df, indices=cells, fields="intensities", drop=TRUE, verbose=verbose);
    verbose && str(verbose, y);
    verbose && exit(verbose);
    # Not needed anymore
    cells <- NULL;

    # Shift?
    if (shift != 0) {
      verbose && enter(verbose, "Shifting signals");
      y <- y + shift;
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Estimating parameters");
    mu <- median(y, na.rm=TRUE);
    # Not needed anymore
    y <- NULL;
    verbose && exit(verbose);

    fitKK <- list(mu=mu, unitType=unitType, nbrOfUnits=nbrOfUnits, nbrOfCells=nbrOfCells);
    verbose && str(verbose, fitKK);

    fit[[key]] <- fitKK;
    # Not needed anymore
    fitKK <- NULL;
    verbose && exit(verbose);
  } # for (kk ...)

  verbose && cat(verbose, "Fitted parameters:");
  verbose && str(verbose, fit);

  verbose && exit(verbose);

  fit;
}, protected=TRUE)





setMethodS3("getNormalizeSignalsOne", "UnitTypeScaleNormalization", function(this, df, fit, ..., verbose=FALSE) {
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
  units <- params$unitsToUpdate;
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  shift <- params$shift;
  stratifyBy <- params$typesToFit;
  targetAvg <- params$targetAvg;
  verbose && exit(verbose);

  # Private/temporary/test code /HB 2008-07-28
  if (isTRUE(this$.foo)) {
    verbose && enter(verbose, "Getting target average for this array");
    unitTypes <- sapply(fit, FUN=function(x) x$unitType);
    mus <- sapply(fit, FUN=function(x) x$mu);
    keep <- (unitTypes %in% c(2,5));
    mus <- mus[keep];
    targetAvg <- mean(mus);
    verbose && cat(verbose, "targetAvg: ", targetAvg);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Getting unit types");
  cdf <- getCdf(df);
  unitTypes <- getUnitTypes(cdf, units=units, verbose=verbose);
  uniqueUnitTypes <- sort(unique(unitTypes));
  knownUnitTypes <- attr(unitTypes, "types");
  verbose && exit(verbose);


  data <- data.frame(cell=NULL, y=NULL);
  for (kk in seq_along(uniqueUnitTypes)) {
    unitType <- uniqueUnitTypes[kk];
    key <- names(knownUnitTypes)[unitType == knownUnitTypes];
    names(unitType) <- key;
    verbose && enter(verbose, sprintf("Unit type %d (%d; '%s') of %d",
                    kk, unitType, key, length(uniqueUnitTypes)));
    unitsKK <- units[unitTypes == unitType];
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, unitsKK);
    verbose && exit(verbose);

    verbose && enter(verbose, "Identifying cell indices");
    cells <- getCellsToInternal(this, units=unitsKK, stratifyBy=stratifyBy, verbose=verbose);
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading signals");
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
    mu <- fit[[key]]$mu;
    b <- targetAvg / mu;
    verbose && printf(verbose, "Scale factor: %.2f\n", b);
    y <- b*y;
    verbose && cat(verbose, "Normalized signals:");
    verbose && str(verbose, y);
    # Sanity check
    yM <- median(y, na.rm=TRUE);
    verbose && printf(verbose, "Median after: %.2f\n", yM);
    verbose && exit(verbose);

    dataKK <- data.frame(cell=cells, y=y);
    # Not needed anymore
    cells <- y <- mu <- b <- yM <- NULL;

    data <- rbind(data, dataKK);
    # Not needed anymore
    dataKK <- NULL;
    verbose && exit(verbose);
  } # for (kk ...)


  verbose && cat(verbose, "Normalized data:");
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
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
setMethodS3("process", "UnitTypeScaleNormalization", function(this, ..., skip=FALSE, force=FALSE, verbose=FALSE) {
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
    data <- getNormalizeSignalsOne(this, df=df, fit=fit,
                                                  verbose=less(verbose));
    verbose && str(verbose, data);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store normalized data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Writing normalized data");
    writeSignals(this, pathname=pathname, cells=data$cell,
             intensities=data$y, templateFile=df, verbose=less(verbose));
    # Not needed anymore
    data <- NULL;

    ## Create checksum file
    dfZ <- getChecksumFile(pathname)

    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

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
# 2008-07-26
# o Created from ScaleNormalization3.R.
############################################################################
