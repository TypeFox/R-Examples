###########################################################################/**
# @RdocClass ScaleNormalization
#
# @title "The ScaleNormalization class"
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
#   \item{dataSet}{@see "AffymetrixCelSet" to be normalized.}
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{targetAvg}{A @numeric value.}
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{subsetToAvg}{The probes to calculate average signal over.
#     If a single @numeric in (0,1), then this
#     fraction of all probes will be used.
#     If @NULL, all probes are considered.}
#   \item{typesToAvg}{Types of probes to be used when calculating the
#     average signal.
#     If \code{"pm"} and \code{"mm"} only perfect-match and mismatch
#     probes are used, respectively. If \code{"pmmm"} both types are used.
#   }
#   \item{shift}{Optional amount of shift if data before fitting/normalizing.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("ScaleNormalization", function(dataSet=NULL, ..., targetAvg=4400, subsetToUpdate=NULL, typesToUpdate=NULL, subsetToAvg="-XY", typesToAvg=typesToUpdate, shift=0) {
  extraTags <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(dataSet)) {
    # Argument 'targetAvg':
    targetAvg <- Arguments$getDouble(targetAvg, range=c(1,Inf));

    # Argument 'subsetToAvg':
    if (is.null(subsetToAvg)) {
    } else if (is.character(subsetToAvg)) {
      if (subsetToAvg %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'subsetToAvg': ", subsetToAvg);
      }
      extraTags <- c(extraTags, subsetToAvg=subsetToAvg);
    } else {
      cdf <- getCdf(dataSet);
      subsetToAvg <- Arguments$getIndices(subsetToAvg, max=nbrOfUnits(cdf));
      subsetToAvg <- unique(subsetToAvg);
      subsetToAvg <- sort(subsetToAvg);
    }

    # Argument 'shift':
    shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"));
  }

  extend(ProbeLevelTransform(dataSet=dataSet, ...), "ScaleNormalization",
    .subsetToUpdate = subsetToUpdate,
    .typesToUpdate = typesToUpdate,
    .targetAvg = targetAvg,
    .subsetToAvg = subsetToAvg,
    .typesToAvg = typesToAvg,
    .extraTags = extraTags,
    shift = shift
  )
})


setMethodS3("getAsteriskTags", "ScaleNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Add class-specific tags
  shift <- as.integer(round(this$shift));
  if (shift != 0) {
    tags <- c(tags, sprintf("%+d", shift));
  }

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getSubsetToUpdate", "ScaleNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  subsetToUpdate <- this$.subsetToUpdate;

  # Done?
  if (identical(attr(subsetToUpdate, "adjusted"), TRUE))
    return(subsetToUpdate);

  # Ad hoc solution for ChipEffectSet:s for now. /HB 2007-04-11
  ds <- getInputDataSet(this);
  if (inherits(ds, "ChipEffectSet")) {
    verbose && enter(verbose, "Identifying possible cells in ", class(ds)[1]);
    df <- getOneFile(ds);
    # Cannot use 'unlist=TRUE' next, because restructuring might occur.
    possibleCells <- getCellIndices(df, verbose=less(verbose));
    possibleCells <- unlist(possibleCells, use.names=FALSE);
    possibleCells <- sort(possibleCells);
    verbose && str(verbose, possibleCells);

    verbose && cat(verbose, "'subsetToUpdate' (before): ");
    verbose && str(verbose, subsetToUpdate);

    if (is.null(subsetToUpdate)) {
      subsetToUpdate <- possibleCells;
    } else {
      subsetToUpdate <- intersect(subsetToUpdate, possibleCells);
    }
    # Not needed anymore
    possibleCells <- NULL;
    verbose && cat(verbose, "'subsetToUpdate' (after): ");
    verbose && str(verbose, subsetToUpdate);

    attr(subsetToUpdate, "adjusted") <- TRUE;
    this$.subsetToUpdate <- subsetToUpdate;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  }

  subsetToUpdate;
}, private=TRUE)




setMethodS3("getSubsetToAvg", "ScaleNormalization", function(this, ..., verbose=FALSE) {
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

  # Done?
  if (identical(attr(subsetToAvg, "adjusted"), TRUE))
    return(subsetToAvg);


  ds <- getInputDataSet(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset with a prespecified set of units?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.character(subsetToAvg)) {
    if (subsetToAvg %in% c("-X", "-Y", "-XY")) {
      verbose && enter(verbose, "Identify subset of units from genome information");
      verbose && cat(verbose, "subsetToAvg: ", subsetToAvg);

      # Look up in cache
      subset <- this$.subsetToAvgExpanded;
      if (is.null(subset)) {
        cdf <- getCdf(ds);
        gi <- getGenomeInformation(cdf);
        # Get the genome information (throws an exception if missing)
        verbose && print(verbose, gi);

        # Identify units to be excluded
        if (subsetToAvg == "-X") {
          units <- getUnitsOnChromosomes(gi, 23, .checkArgs=FALSE);
        } else if (subsetToAvg == "-Y") {
          units <- getUnitsOnChromosomes(gi, 24, .checkArgs=FALSE);
        } else if (subsetToAvg == "-XY") {
          units <- getUnitsOnChromosomes(gi, 23:24, .checkArgs=FALSE);
        }

        verbose && cat(verbose, "Units to exclude: ");
        verbose && str(verbose, units);

        # The units to keep
        units <- setdiff(1:nbrOfUnits(cdf), units);

        verbose && cat(verbose, "Units to include: ");
        verbose && str(verbose, units);

        # Identify cell indices for these units
        subset <- getCellIndices(cdf, units=units, unlist=TRUE);

        # Store
        this$.subsetToAvgExpanded <- subset;
      } # if (is.null(subset))

      subsetToAvg <- subset;
      # Not needed anymore
      subset <- NULL;

      verbose && exit(verbose);
    }
  }


  # Ad hoc solution for ChipEffectSet:s for now. /HB 2007-04-11
  ds <- getInputDataSet(this);
  if (inherits(ds, "ChipEffectSet")) {
    verbose && enter(verbose, "Identifying possible cells in ", class(ds)[1]);
    df <- getOneFile(ds);
    # Cannot use 'unlist=TRUE' next, because restructuring might occur.
    possibleCells <- getCellIndices(df, verbose=less(verbose));
    possibleCells <- unlist(possibleCells, use.names=FALSE);
    possibleCells <- sort(possibleCells);
    verbose && str(verbose, possibleCells);

    verbose && cat(verbose, "'subsetToAvg' (before): ");
    verbose && str(verbose, subsetToAvg);

    if (is.null(subsetToAvg)) {
      subsetToAvg <- possibleCells;
    } else {
      subsetToAvg <- intersect(subsetToAvg, possibleCells);
    }
    # Not needed anymore
    possibleCells <- NULL;
    verbose && cat(verbose, "'subsetToAvg' (after): ");
    verbose && str(verbose, subsetToAvg);

    attr(subsetToAvg, "adjusted") <- TRUE;
    this$.subsetToAvg <- subsetToAvg;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  }

  subsetToAvg;
}, private=TRUE)




setMethodS3("getParameters", "ScaleNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    subsetToUpdate = getSubsetToUpdate(this),
    typesToUpdate = this$.typesToUpdate,
    subsetToAvg = getSubsetToAvg(this),
    typesToAvg = this$.typesToAvg,
    targetAvg = this$.targetAvg,
    shift = this$shift
  );

  # Append the two sets
  params <- c(params, params2);

  params;
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
setMethodS3("process", "ScaleNormalization", function(this, ..., skip=FALSE, force=FALSE, verbose=FALSE) {
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
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);

  # Get algorithm parameters
  params <- getParameters(this);
  targetAvg <- params$targetAvg;

  subsetToAvg <- getSubsetToAvg(this, verbose=less(verbose));
  verbose && cat(verbose, "subsetToAvg:");
  verbose && str(verbose, subsetToAvg);
  verbose && print(verbose, range(subsetToAvg));

  subsetToUpdate <- getSubsetToUpdate(this, verbose=less(verbose));
  verbose && cat(verbose, "subsetToUpdate:");
  verbose && str(verbose, subsetToUpdate);
  verbose && print(verbose, range(subsetToUpdate));

  # Get the output path
  outputPath <- getPath(this);


  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing ", length(ds), " arrays");
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
      verbose && exit(verbose);
      next;
    }

    # Get probe signals for the fit
    verbose && enter(verbose, "Reading probe intensities for fitting");
##    data <- getDataFlat(ce, units=map, fields="theta", verbose=less(verbose));
    x <- getData(df, fields="intensities", indices=subsetToAvg,
                                 verbose=less(verbose,2))$intensities;
    verbose && str(verbose, x);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Shift?
    shift <- this$shift;
    if (shift != 0) {
      x <- x + shift;
    }

    # Estimating scale
    verbose && enter(verbose, "Estimating scale");
    xM <- median(x, na.rm=TRUE);
    verbose && printf(verbose, "Median before: %.2f\n", xM);
    verbose && printf(verbose, "Target avg: %.2f\n", targetAvg);
    b <- targetAvg / xM;
    verbose && printf(verbose, "Scale: %.2f\n", b);
    verbose && exit(verbose);

    # Garbage collect
    # Not needed anymore
    x <- xM <- NULL;
    gc <- gc();
    verbose && print(verbose, gc);

    # Get probe signals to be updated
    verbose && enter(verbose, "Getting signals");
##    data <- getDataFlat(ce, units=map, fields="theta", verbose=less(verbose));
    x <- getData(df, fields="intensities", indices=subsetToUpdate,
                                    verbose=less(verbose,2))$intensities;
    verbose && str(verbose, x);
    verbose && exit(verbose);

    # Rescale
    verbose && enter(verbose, "Rescaling");
    x <- b*x;
    verbose && str(verbose, x);
    xM <- median(x, na.rm=TRUE);
    verbose && printf(verbose, "Median after: %.2f\n", xM);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Write normalized data to file
    verbose && enter(verbose, "Writing normalized probe signals");

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing");

    # Write to a temporary file (allow rename of existing one if forced)
    isFile <- (!skip && isFile(pathname));
    pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

    createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose));
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing normalized signals");
#    updateDataFlat(ceN, data=data, verbose=less(verbose));
#    # Not needed anymore
#    data <- NULL;
    .updateCel(pathnameT, indices=subsetToUpdate, intensities=x);
    # Not needed anymore
    x <- NULL;

    # Rename temporary file
    popTemporaryFile(pathnameT, verbose=verbose);

    ## Create checksum file
    dfZ <- getChecksumFile(pathname)

    verbose && exit(verbose);

    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  # Create result set
  outputDataSet <- getOutputDataSet(this, force=TRUE, verbose=less(verbose));

  # Update the output data set
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);

  outputDataSet;
})

############################################################################
# HISTORY:
# 2008-02-19
# o Added support for "-X", "-Y" and "-XY" for argument 'subsetToAvg'.
# o Added support for constructor argument 'shift'.
# 2007-09-18
# o Updated all getCellIndices() to use 'unlist=TRUE'.
# 2007-04-16
# o Created.
############################################################################
