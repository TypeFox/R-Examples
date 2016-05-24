###########################################################################/**
# @set "class=ProbeLevelModel"
# @RdocMethod fit
#
# @title "Estimates the model parameters"
#
# \description{
#  @get "title" for all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be fitted.
#     If @NULL, all units are considered.
#     If \code{remaining}, only non-fitted units are considered.
#   }
#   \item{...}{Arguments passed to \code{readUnits()}.}
#   \item{force}{If @TRUE, already fitted units are re-fitted, and
#     cached data is re-read.}
#   \item{ram}{A @double indicating if more or less units should
#     be loaded into memory at the same time.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @integer @vector of indices of the units fitted,
#  or @NULL if no units was (had to be) fitted.
# }
#
# \details{
#  All estimates are stored to file.
#
#  The non-array specific parameter estimates together with standard deviation
#  estimates and convergence information are stored in one file.
#
#  The parameter estimates specific to each array, typically "chip effects",
#  are stored in array specific files.
#
#   Data set specific estimates [L = number of probes]:
#    phi [L doubles] (probe affinities), sd(phi) [L doubles],
#    isOutlier(phi) [L logicals]
#
#   Algorithm-specific results:
#    iter [1 integer], convergence1 [1 logical], convergence2 [1 logical]
#    dTheta [1 double]
#    sd(eps) - [1 double] estimated standard deviation of the error term
#
#   Array-specific estimates [K = nbr of arrays]:
#    theta [K doubles] (chip effects), sd(theta) [K doubles],
#    isOutlier(theta) [K logicals]
#
#   For each array and each unit group, we store:
#     1 theta, 1 sd(theta), 1 isOutlier(theta), i.e. (float, float, bit)
#   => For each array and each unit (with \eqn{G_j} groups), we store:
#     \eqn{G_j} theta, \eqn{G_j} sd(theta), \eqn{G_j} isOutlier(theta),
#   i.e. \eqn{G_j}*(float, float, bit).
#   For optimal access we store all thetas first, then all sd(theta) and the
#   all isOutlier(theta).
#   To keep track of the number of groups in each unit, we have to have a
#   (unit, ngroups) map.  This can be obtained from getUnitNames() for the
#   AffymetrixCdfFile class.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fit", "ProbeLevelModel", function(this, units="remaining", ..., force=FALSE, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the some basic information about this model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cs <- getDataSet(this);
  cdf <- getCdf(cs);
  nbrOfArrays <- length(cs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  doRemaining <- FALSE;
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  } else if (identical(units, "remaining")) {
    doRemaining <- TRUE;
  } else {
    throw("Unknown mode of argument 'units': ", mode(units));
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Fitting model of class ", class(this)[1]);

  verbose && print(verbose, this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify units to be fitted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    units <- 1:nbrOfUnits;
  } else if (doRemaining) {
    verbose && enter(verbose, "Identifying non-estimated units")
    units <- findUnitsTodo(this, verbose=less(verbose));
    nbrOfUnits <- length(units);
    verbose && exit(verbose);
  } else {
    # Fit only unique units
    units <- unique(units);
    nbrOfUnits <- length(units);
  }
  verbose && printf(verbose, "Getting model fit for %d units.\n", nbrOfUnits);

  # Identify which of the requested units have *not* already been estimated
  if (!doRemaining) {
    if (force) {
      verbose && printf(verbose, "All of these are forced to be fitted.\n");
    } else {
      units <- findUnitsTodo(this, units=units, verbose=less(verbose));
      nbrOfUnits <- length(units);
      verbose && printf(verbose, "Out of these, %d units need to be fitted.\n", nbrOfUnits);
    }
  }

  # Nothing more to do?
  if (nbrOfUnits == 0)
    return(invisible(NULL));



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for prior parameter estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  priors <- getListOfPriors(this, verbose=log);
  hasPriors <- (length(priors) > 0);
  if (hasPriors) {
    verbose && cat(verbose, "Prior parameters detected");
  }
  # Not needed anymore
  priors <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit all units that with single-cell unit groups using special
  # fit function, if available. This can speed up the fitting
  # substantially.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##   fitSingleCellUnit <- getFitSingleCellUnitFunction(this);
##   if (!is.null(fitSingleCellUnit)) {
##     verbose && enter(verbose, "Fitting single-cell units");
##
##     verbose && enter(verbose, "Identifying single-cell units");
##     counts <- nbrOfCellsPerUnit(cdf, units=units, verbose=less(verbose, 5));
##     verbose && print(verbose, table(verbose));
##     singleCellUnits <- which(counts == 1);
##     # Not needed anymore
##     counts <- NULL;
##     verbose && str(verbose, singleCellUnits);
##     verbose && exit(verbose);
##
##     # Nothing to do?
##     if (length(singleCellUnits) > 0) {
##       verbose && enter(verbose, "Identifying CDF cell indices");
##       cells <- getCellIndices(this, units=units, unlist=TRUE, useNames=FALSE, ...);
##       verbose && str(verbose, cells);
##       # Sanity check
##       stopifnot(length(cells) == length(singleCellUnits));
##       verbose && exit(verbose);
##       verbose && exit(verbose);
##
##       verbose && enter(verbose, "Reading signals");
##       verbose && exit(verbose);
##       # Not needed anymore
##       cells <- NULL;
##
##       verbose && enter(verbose, "Fitting units");
##       verbose && exit(verbose);
##
##       verbose && enter(verbose, "Storing estimates");
##       verbose && exit(verbose);
##
##       # Update remaining units to do
## ##      units <- setdiff(units, singleCellUnits);
##       nbrOfUnits <- length(units);
##     }
##
##     verbose && exit(verbose);
##
##     # Nothing more to do?
##     if (nbrOfUnits == 0)
##       return(NULL);
##   }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify unit types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying unit types:");
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  unitTypes <- getUnitTypes(cdf, units=units);
  verbose && cat(verbose, "Unit types:");
  verbose && str(verbose, unitTypes);
  verbose && print(verbose, table(unitTypes));
  uUnitTypes <- sort(unique(unitTypes));
  verbose && cat(verbose, "Unique unit types:");
  types <- attr(unitTypes, "types");
  names(uUnitTypes) <- names(types)[match(uUnitTypes, types)];
  verbose && print(verbose, uUnitTypes);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting parameter sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up parameter sets");
  # Get the model-fit function
  fitUnit <- getFitUnitFunction(this);

  # Get (and create if missing) the probe-affinity file (one per data set)
  paf <- getProbeAffinityFile(this, verbose=less(verbose));

  # Get (and create if missing) the chip-effect files (one per array)
  ces <- getChipEffectSet(this, ram=ram, verbose=less(verbose));

  # Garbage collect
  clearCache(cs);
  clearCache(paf);
  clearCache(ces);
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit each type of units independently
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting ", class(this)[1], " for each unit type separately");
  verbose && cat(verbose, "Unit types:");
  verbose && print(verbose, uUnitTypes);

  for (tt in seq_along(uUnitTypes)) {
    unitType <- uUnitTypes[tt];
    unitTypeLabel <- names(uUnitTypes)[tt];
    verbose && enter(verbose, sprintf("Unit type #%d ('%s') of %d", tt, unitTypeLabel, length(uUnitTypes)));

    unitsTT <- units[unitTypes == unitType];
    nbrOfUnitsTT <- length(unitsTT);
    verbose && printf(verbose, "Unit type: %s (code=%d)\n", unitTypeLabel, unitType);
    verbose && cat(verbose, "Number of units of this type: ", nbrOfUnitsTT);
    verbose && str(verbose, unitsTT);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit the model by unit dimensions (at least for the large classes)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Fitting the model by unit dimensions (at least for the large classes)");

    verbose && enter(verbose, "Grouping units into equivalent (unit,group,cell) dimensions");
    unitDimensions <- groupUnitsByDimension(cdf, units=unitsTT, verbose=less(verbose, 50));
    # Not needed anymore
    # Not needed anymore
    unitsTT <- NULL;

    sets <- unitDimensions$sets;
    dims <- unitDimensions$setDimensions;
    o <- order(dims$nbrOfUnits, decreasing=TRUE);
    dims <- dims[o,];
    sets <- sets[o];

    if (verbose) {
      verbose && printf(verbose, "%d classes of unit dimensions:\n", nrow(dims));
      dimsT <- dims;
      dimsT$nbrOfCellsPerGroup <- sapply(dims$nbrOfCellsPerGroup, FUN=hpaste, collapse="+");
      if (nrow(dimsT) < 100) {
        verbose && print(verbose, dimsT);
      } else {
        verbose && print(verbose, head(dimsT));
        verbose && print(verbose, tail(dimsT));
      }
      # Not needed anymore
      dimsT <- NULL;
    }
    verbose && exit(verbose);


    # Identify large classes
    minNbrOfUnits <- getOption(aromaSettings, "models/Plm/chunks/minNbrOfUnits", 500L);
    isLarge <- (dims$nbrOfUnits >= minNbrOfUnits);
    nLarge <- sum(isLarge);
    nTotal <- nrow(dims);
    nuLarge <- sum(dims$nbrOfUnits[isLarge]);
    nuTotal <- sum(dims$nbrOfUnits);
    verbose && printf(verbose, "There are %d large classes of units, each with a unique dimension and that each contains at least %d units.  In total these large classes constitute %d (%.2f%%) units.\n", nLarge, minNbrOfUnits, nuLarge, 100*nuLarge/nuTotal);

    large <- mapply(sets[isLarge], dims$nbrOfCellsPerGroup[isLarge], FUN=function(set, dim) {
      list(units=set$units, dim=dim);
    }, SIMPLIFY=FALSE);
    names(large) <- sapply(dims$nbrOfCellsPerGroup[isLarge], FUN=paste, collapse="+");
    dimChunks <- large;

    # Additional small sets, if any
    small <- lapply(sets[!isLarge], FUN=function(x) x$units);
    small <- unlist(small, use.names=FALSE);
    if (length(small) > 0) {
      small <- list(mix=list(units=small));
      dimChunks <- c(dimChunks, small);
    }

    for (kk in seq_along(dimChunks)) {
      dimChunk <- dimChunks[[kk]];
      dim <- dimChunk$dim;
      if (is.null(dim)) {
        dimDescription <- "various dimensions";
      } else {
        nbrOfGroups <- length(dim);
        dimStr <- paste(dimChunk$dim, collapse="+");
        dimDescription <- sprintf("%d groups/%s cells", nbrOfGroups, dimStr);
      }

      verbose && enter(verbose, sprintf("Unit dimension #%d (%s) of %d", kk, dimDescription, length(dimChunks)));

      unitsKK <- dimChunk$units;
      nbrOfUnitsKK <- length(unitsKK);
      verbose && printf(verbose, "Number of units: %d (%.2f%%) out of %d '%s' units (code=%d)\n", nbrOfUnitsKK, 100*nbrOfUnitsKK/nbrOfUnitsTT, nbrOfUnitsTT, unitTypeLabel, unitType);
      verbose && str(verbose, unitsKK);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit the model in chunks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Calculating number of units to fit per chunk");

    verbose && cat(verbose, "RAM scale factor: ", ram);

    bytesPerChunk <- 100e6;       # 100Mb
    verbose && cat(verbose, "Bytes per chunk: ", bytesPerChunk);

    bytesPerUnitAndArray <- 500;  # Just a rough number; good enough?
    verbose && cat(verbose, "Bytes per unit and array: ", bytesPerUnitAndArray);

    bytesPerUnit <- nbrOfArrays * bytesPerUnitAndArray;
    verbose && cat(verbose, "Bytes per unit: ", bytesPerUnit);

    unitsPerChunk <- ram * bytesPerChunk / bytesPerUnit;
    unitsPerChunk <- as.integer(max(unitsPerChunk,1));
    verbose && cat(verbose, "Number of units per chunk: ", unitsPerChunk);

    nbrOfChunks <- ceiling(nbrOfUnitsKK / unitsPerChunk);
    verbose && cat(verbose, "Number of chunks: ", nbrOfChunks);

    idxs <- 1:nbrOfUnitsKK;
    head <- 1:unitsPerChunk;

    verbose && exit(verbose);


    # Time the fitting.
    startTime <- processTime();

    timers <- list(total=0, read=0, fit=0, writePaf=0, writeCes=0, gc=0);

    count <- 1;
    while (length(idxs) > 0) {
      tTotal <- processTime();

      verbose && enter(verbose, sprintf("Fitting chunk #%d of %d of '%s' units (code=%d) with %s", count, nbrOfChunks, unitTypeLabel, unitType, dimDescription));
      if (length(idxs) < unitsPerChunk) {
        head <- 1:length(idxs);
      }
      uu <- idxs[head];

      unitsChunk <- unitsKK[uu];
      verbose && printf(verbose, "Unit type: '%s' (code=%d)\n", unitTypeLabel, unitType);
      verbose && cat(verbose, "Unit dimension: ", dimDescription);
      verbose && cat(verbose, "Number of units in chunk: ", length(unitsChunk));
      verbose && cat(verbose, "Units: ");
      verbose && str(verbose, unitsChunk);

      # Sanitcy check
      stopifnot(length(idxs) > 0);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Get the CEL intensities by units
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tRead <- processTime();
      y <- readUnits(this, units=unitsChunk, ..., force=force, cache=FALSE,
                                                      verbose=less(verbose));
      timers$read <- timers$read + (processTime() - tRead);

##    # It is not possible to do the below sanity check, because
##    # readUnits() may have merged groups etc in a way we cannot
##    # know.  /HB 2012-01-11
##      # Sanity checks
##      if (!is.null(dim)) {
##        nbrOfGroups <- length(dim);
##        if (length(y[[1]]) != nbrOfGroups) {
##          throw("Internal error: The read units does not have ", nbrOfGroups, " groups: ", length(y[[1]]));
##        }
##      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Read prior parameter estimates
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (hasPriors) {
        priors <- readPriorsByUnits(this, units=unitsChunk, ..., force=force,
                                         cache=FALSE, verbose=less(verbose));
        timers$read <- timers$read + (processTime() - tRead);
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Fit the model
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Fitting probe-level model");
      tFit <- processTime();

      if (hasPriors) {
        verbose && cat(verbose, "Calling fitUnit() via mapply() with prior parameter estimates");
        fit <- mapply(y, priors=priors, FUN=fitUnit, SIMPLIFY=FALSE);
      } else {
        verbose && cat(verbose, "Calling fitUnit() via lapply()");
        fit <- lapply(y, FUN=fitUnit);
      }

      timers$fit <- timers$fit + (processTime() - tFit);
      y <- NULL; # Not needed anymore (to minimize memory usage)
      verbose && str(verbose, fit[1]);
      verbose && exit(verbose);

      # Garbage collection
      tGc <- processTime();
      gc <- gc();
      timers$gc <- timers$gc + (processTime() - tGc);
      verbose && print(verbose, gc);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store probe affinities
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing probe-affinity estimates");
      tWritePaf <- processTime();
      updateUnits(paf, units=unitsChunk, data=fit, verbose=less(verbose));
      timers$writePaf <- timers$writePaf + (processTime() - tWritePaf);
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store chip effects
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing chip-effect estimates");
      tWriteCes <- processTime();
      updateUnits(ces, units=unitsChunk, data=fit, verbose=less(verbose));
      timers$writeCes <- timers$writeCes + (processTime() - tWriteCes);
      verbose && exit(verbose);

      fit <- NULL; # Not needed anymore

      # Next chunk
      idxs <- idxs[-head];
      count <- count + 1;

      # Garbage collection
      tGc <- processTime();
      gc <- gc();
      verbose && print(verbose, gc);
      timers$gc <- timers$gc + (processTime() - tGc);

      timers$total <- timers$total + (processTime() - tTotal);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # ETA
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (verbose) {
        # Clarifies itself once in a while (in case running two in parallel).
        verbose && print(verbose, this);

        # Fraction left
        fLeft <- length(idxs) / nbrOfUnits;
        # Time this far
        lapTime <- processTime() - startTime;
        t <- Sys.time() - lapTime[3];
        printf(verbose, "Started: %s\n", format(t, "%Y%m%d %H:%M:%S"));
        # Estimated time left
        fDone <- 1-fLeft;
        timeLeft <- fLeft/fDone * lapTime;
        t <- timeLeft[3];
        printf(verbose, "Estimated time left for unit type '%s': %.1fmin\n", unitTypeLabel, t/60);
        # Estimate time to arrivale
        eta <- Sys.time() + t;
        printf(verbose, "ETA for unit type '%s': %s\n", unitTypeLabel, format(eta, "%Y%m%d %H:%M:%S"));
      }

      verbose && exit(verbose);
    } # while(length(idxs) > 0)  # For each chunk of units...

    totalTime <- processTime() - startTime;
    if (verbose) {
      nunits <- length(unitsKK);
      t <- totalTime[3];
      printf(verbose, "Total time for all '%s' units across all %d arrays: %.2fs == %.2fmin\n", unitTypeLabel, nbrOfArrays, t, t/60);
      t <- totalTime[3]/nunits;
      printf(verbose, "Total time per '%s' unit across all %d arrays: %.2fs/unit\n", unitTypeLabel, nbrOfArrays, t);
      t <- totalTime[3]/nunits/nbrOfArrays;
      printf(verbose, "Total time per '%s' unit and array: %.3gms/unit & array\n", unitTypeLabel, 1000*t);
      t <- nbrOfUnits(cdf)*totalTime[3]/nunits/nbrOfArrays;
      printf(verbose, "Total time for one array (%d units): %.2fmin = %.2fh\n", nbrOfUnits(cdf), t/60, t/3600);
      t <- nbrOfUnits(cdf)*totalTime[3]/nunits;
      printf(verbose, "Total time for complete data set: %.2fmin = %.2fh\n", t/60, t/3600);
      # Get distribution of what is spend where
      timers$write <- timers$writePaf + timers$writeCes;
      t <- lapply(timers, FUN=function(timer) unname(timer[3]));
      t <- unlist(t);
      t <- 100 * t / t["total"];
      printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%% (of which %.2f%% is for encoding/writing chip-effects), Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], 100*t["writeCes"]/t["write"], t["gc"]);
    }

      verbose && exit(verbose);
    } # for (kk ...) # For each unit dimension class...
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (tt ...)  # For each unit types...


  ## Generate checksum files
  pafZ <- getChecksumFile(paf)
  cesZ <- getChecksumFileSet(ces)

  verbose && exit(verbose);

  invisible(units);
}) # fit()




############################################################################
# HISTORY:
# 2012-01-11
# o BUG FIX: fit() for ProbeLevelModel would throw an error due to a
#   sanity check on unit dimensions that was only valid in certain cases.
# 2011-11-20
# o BUG FIX: In the most recent update to fit() when fitting "remaining"
#   sets of units ("various dimensions") would refit everything iff there
#   where no such units (leading to units == NULL).
# 2011-11-18
# o Now fit() for ProbeLevelModel fits one type of units at the time,
#   which in turn is fitted in chunks of units with equal number of
#   groups and cells per groups (very small chunks are merged together).
# 2011-11-10
# o Now fit() for ProbeLevelModel fits one type of units at the time.
#   For each unit type, the fitting is done in chunks.  The timing
#   statistics presented in the verbose output is now per unit type.
# 2010-02-16
# o Added additional startup verbose output in fit() for ProbeLevelModel.R.
# 2009-02-21
# o Removed deprecated argument 'moreUnits' of fit() of ProbeLevelModel.
# 2008-12-08
# o Added basic support for priors in fit() of ProbeLevelModel.
# 2008-07-22
# o Now argument 'ram' is passed down to getChipEffectSet() which in turn
#   pass it down to getMonocellCdf(), which pass it to createMonocellCdf()
#   in case the monocell CDF is missing.  This will increase the chances
#   that fit(..., ram=<small value>) will work with small amount of RAM.
# 2008-05-31
# o Removed an obsolete debug print() statement.
# 2008-02-17
# o Moved fit() to its own source file.
# 2007-12-10
# o Now getChipEffectSet() of ProbeLevelModel infers the monocell CDF from
#   the CDF of the input data set and uses that when retrieving the
#   chip-effect CEL set.  In other words, if the CDF is overridden for
#   the input data set, it will also be overridden (with the corresponding
#   monocell CDF) in the chip-effect set.  Before the monocell CDF was
#   always inferred from the CEL header, if the CEL file existed.
# 2007-12-08
# o Now the tag for the 'shift' is also set in getAsteriskTag().
# o Now the tag for the 'probeModel' is set in getAsteriskTag().
# o Added argument 'shift' to ProbeLevelModel.
# 2007-08-09
# o getProbeAffinityFile() of ProbeLevelModel now creates CEL files with
#   upper-case filename extension "*.CEL", not "*.cel".  The reason for this
#   is that some software don't recognize lower case filename extensions :(
# 2007-04-12
# o Added 'force' argument to getResidualSet().
# 2007-04-02
# o Added support for the "pm+mm" probe model.
# 2007-02-29
# o BUG FIX: Probe-affinities was not save, resulting in all zeroes.
#   This was due to renaming getProbeAffinites() to getProbeAffinityFile().
# 2007-02-28
# o Added ETA to verbose output of fit() for the ProbeLevelModel.
# o Memory optimization: Further memory optimization by clearing the
#   cache of the 'cs', the 'paf', and the 'ces', before fitting.
# 2007-02-22
# o Added getChipEffectSet() and getProbeAffinityFile() to replace
#   getChipEffects() and getProbeAffinites() in some future version.
# 2007-02-09
# o Added an additional garbage collection after fitting the PLM, but
#   before storing parameter estimates.
# 2007-01-06
# o Now gc() memory information is outputted after each chunk.
# o Updated the formula for calculating the number of units per chunk in
#   fit(). Gives roughly the same number of units though.
# o Added probe model 'min1(PM-MM)' for modelling y = min(PM-MM,1), which
#   is how dChip v2006-12-14 is doing it.
# o Now ProbeLevelModel inherits directly from UnitModel.
# 2007-01-03
# o "Protected" several methods to simplify the interface for end users.
# o Added support from "PM-MM" probe models in addition to the default
#   "PM only" model.
# 2006-09-26
# o Fixed the timers for fit(). They only worked so and so before and only
#   only Windows.  Using processTime()
# 2006-09-14
# o Added detailed timing information to the verbose output of fit().
# 2006-09-11
# o Added argument ProbeLevelModel(..., standardize=TRUE) to make the
#   results from different PLMs be on the same scale.
# 2006-09-10
# o Added findUnitsTodo().
# o Update fit() to make use of new ChipEffects and ProbeAffinity classes.
#   The code is much cleaner now!
# 2006-08-28
# o Added plotMvsPosition().
# 2006-08-26
# o Added new getChipEffects().
# 2006-08-25
# o Created from AffymetrixLiWongModel.  So much is in common with the RMA
#   model so a lot can be reused if done cleverly.
############################################################################
