setMethodS3("boxplotStats", "ChipEffectSet", function(this, type=c("NUSE", "RLE"), transform=NULL, ...) {
  if (toupper(type) == "NUSE") {
    calculateNuseBoxplotStats(this, ...);
  } else if (toupper(type) == "RLE") {
    calculateRleBoxplotStats(this, ...);
  } else {
    calculateFieldBoxplotStats(this, field=type, transform=transform, ...);
  }
})




setMethodS3("calculateFieldBoxplotStats", "ChipEffectSet", function(this, field=c("theta", "sdTheta"), transform=NULL, arrays=NULL, subset=NULL, ..., merge=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # So that one can call plotRle(qa, show.names=FALSE), i.e. not
  # worry about what arguments are in '...'.
  boxplotStats <- appendVarArgs(boxplot.stats);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(this);
  cdfMono <- getCdf(this);
  nbrOfUnits <- nbrOfUnits(cdfMono);

  # Argument 'field':
  field <- match.arg(field);

  # Argument 'transform':
  if (is.null(transform)) {
  } else if (is.function(transform)) {
  } else {
    throw("Argument 'transform' is not a function: ", class(transform)[1]);
  }


  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq_len(nbrOfArrays);
  } else {
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays);
    nbrOfArrays <- length(arrays);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'subset':
  if (!(is.null(subset))) {
    getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1);
    if (!getFraction) {
      units <- Arguments$getIndices(subset, max=nbrOfUnits);
    } else {
      units <- identifyCells(cdfMono, indices=subset, verbose=less(verbose));
    }
  } else {
    units <- 1:nbrOfUnits;
  }


  # Get the (unit, group, cell) map
  ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose, 5));

  verbose && enter(verbose, "Calculating '", field,
                  "' statistics for ", nbrOfArrays, " (specified) arrays");
  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);

  # For each file, calculate boxplot statistics
  stats <- list();
  for (kk in seq_along(arrays)) {
    array <- arrays[kk];
    cef <- this[[array]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                          kk, getName(cef), nbrOfArrays));
    data <- extractMatrix(cef, fields=field, units=ugcMap);
    data <- as.vector(data);
    if (is.function(transform)) {
      data <- transform(data);
    }
    stats[[kk]] <- boxplotStats(data, ...);
    # Not needed anymore
    data <- NULL;
    verbose && exit(verbose);
  }
  names(stats) <- getNames(this)[arrays];
  verbose && exit(verbose);

  # Merge boxplot stats?
  if (merge)
    stats <- mergeBoxplotStats(stats);

  attr(stats, "type") <- field;
  attr(stats, "transform") <- transform;

  stats;
}, protected=TRUE) # calculateFieldBoxplotStats()






# ... : additional arguments to bxp().
setMethodS3("calculateRleBoxplotStats", "ChipEffectSet", function(this, arrays=NULL, subset=NULL, ..., merge=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # So that one can call plotRle(qa, show.names=FALSE), i.e. not
  # worry about what arguments are in '...'.
  boxplotStats <- appendVarArgs(boxplot.stats);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(this);
  cdfMono <- getCdf(this);
  nbrOfUnits <- nbrOfUnits(cdfMono);

  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq_len(nbrOfArrays);
  } else {
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays);
    nbrOfArrays <- length(arrays);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'subset':
  if (!(is.null(subset))) {
    getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1);
    if (!getFraction) {
      units <- Arguments$getIndices(subset, max=nbrOfUnits);
    } else {
      units <- identifyCells(cdfMono, indices=subset, verbose=less(verbose));
    }
  } else {
    units <- 1:nbrOfUnits;
  }


  # Get the (unit, group, cell) map
  ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose, 5));

  # get the vector of median stdvs
  verbose && enter(verbose, "Calculating average log chip effects");
  # Calculating the average on the log (but stored on the intensity) scale.
  avg <- getAverageLog(this, field="intensities", verbose=verbose);
  verbose && exit(verbose);

  medianLE <- extractMatrix(avg, field="theta", units=ugcMap,
                                                 verbose=less(verbose, 5));
  medianLE <- log2(as.vector(medianLE));
  # Not needed anymore
  avg <- NULL;

  verbose && enter(verbose, "Calculating RLE statistics for ", nbrOfArrays,
                                                    " (specified) arrays");
  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);

  # For each file, calculate boxplot statistics
  stats <- list();
  for (kk in seq_along(arrays)) {
    array <- arrays[kk];
    cef <- this[[array]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                          kk, getName(cef), nbrOfArrays));

    theta <- extractMatrix(cef, field="theta", units=ugcMap,
                                                verbose=less(verbose, 5));
    theta <- log2(as.vector(theta));

    # Sanity check
    if (length(theta) != length(medianLE)) {
      throw("Internal error: The number of 'theta' does not match the number of 'medianLE': ", length(theta), " != ", length(medianLE));
    }

    stats[[kk]] <- boxplotStats(theta-medianLE, ...);
    # Not needed anymore
    theta <- NULL;
    verbose && exit(verbose);
  }
  # Not needed anymore
  medianLE <- units <- NULL;
  names(stats) <- getNames(this)[arrays];
  verbose && exit(verbose);

  # Merge boxplot stats?
  if (merge)
    stats <- mergeBoxplotStats(stats);

  attr(stats, "type") <- "RLE";

  stats;
}, protected=TRUE) # calculateRleBoxplotStats()



setMethodS3("calculateNuseBoxplotStats", "ChipEffectSet", function(this, arrays=NULL, subset=NULL, ..., merge=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdfMono <- getCdf(this);
  nbrOfUnits <- nbrOfUnits(cdfMono);
  nbrOfArrays <- length(this);

  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq_len(nbrOfArrays);
  } else {
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays);
    nbrOfArrays <- length(arrays);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'subset':
  if (!(is.null(subset))) {
    getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1);
    if (!getFraction) {
      units <- Arguments$getIndices(subset, max=nbrOfUnits);
    } else {
      units <- identifyCells(cdfMono, indices=subset, verbose=less(verbose));
    }
  } else {
    units <- 1:nbrOfUnits;
  }


  # Get the (unit, group, cell) map
  ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose, 5));

  # get the vector of median stdvs
  verbose && enter(verbose, "Extracting average standard errors across all arrays in the set");
  # Calculating the average on the log (but stored on the intensity) scale.
  avg <- getAverageLog(this, field="stdvs", verbose=verbose);
  verbose && exit(verbose);

  # Note, the average of the 'stdvs' is stored in the 'theta' field.
  medianSE <- extractMatrix(avg, field="theta", units=ugcMap,
                                                 verbose=less(verbose, 5));
  # Not needed anymore
  avg <- NULL;
  medianSE <- log2(as.vector(medianSE));

  verbose && enter(verbose, "Calculating NUSE statistics for ", nbrOfArrays,
                                                    " (specified) arrays");
  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);

  # For each file, calculate boxplot statistics
  stats <- list();
  for (kk in seq_along(arrays)) {
    array <- arrays[kk];
    cef <- this[[array]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                          kk, getName(cef), nbrOfArrays));
    stdvs <- extractMatrix(cef, field="sdTheta", units=ugcMap,
                                                verbose=less(verbose, 5));
    stdvs <- log2(as.vector(stdvs));

    # Sanity check
    if (length(stdvs) != length(medianSE)) {
      throw("Internal error: The number of 'stdvs' does not match the number of 'medianSE': ", length(stdvs), " != ", length(medianSE));
    }

    stats[[kk]] <- boxplot.stats(stdvs/medianSE, ...);
    # Not needed anymore
    stdvs <- NULL;
    verbose && exit(verbose);
  }
  # Not needed anymore
  medianSE <- units <- NULL;
  names(stats) <- getNames(this)[arrays];
  verbose && exit(verbose);

  # Merge boxplot stats?
  if (merge)
    stats <- mergeBoxplotStats(stats);

  attr(stats, "type") <- "NUSE";

  stats;
}, protected=TRUE) # calculateNuseStats()


##########################################################################
# HISTORY:
# 2012-08-30 [HB]
# o Now calculateFieldBoxplotStats() and calculateRleBoxplotStats()
#   calls boxplot.stats() via a local boxplotStats() functions that
#   automatically drops unused arguments.
# 2008-02-28 [EP]
# o Now '...' are passed to boxplot.stats().
# o BUG FIX: Now extractMatrix() is used internally to calculate the
#   'RLE' and the 'NUSE' stats.  Before getData() was used which is not
#   "sensitive" to arguments such as 'mergeGroups', 'mergeAlleles' etc.
# 2008-02-25
# o Renamed to make it explicit that it is boxplot stats that are
#   calculated.
# o Now calculate{Nuse|Rle}Stats() returns the list of boxplot stats,
#   not the combined ones.
# o Added boxplotStats() and support calculate{Nuse|Rle}BoxplotStats()
#   adopted from EPs code.
# 2008-02-22
# o Generalized from the QualityAssessmentModel.
##########################################################################
