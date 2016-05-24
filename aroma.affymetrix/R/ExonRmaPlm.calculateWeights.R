
setMethodS3("calculateWeights", "ExonRmaPlm", function(this, units=NULL, ram=NULL, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Lookup MASS::psi.huber() once; '::' is expensive
  MASS_psi.huber <- MASS::psi.huber;

  resFcn <- function(unit, mergeGroups) {
    nbrOfGroups <- length(unit);
    if (mergeGroups) {
      y <- do.call(rbind, lapply(unit, FUN=.subset2, "eps"));
      y <- log2(y);
      madMerged <- 1.4826 * median(abs(y));
    }
    res <- lapply(1:nbrOfGroups, FUN=function(gg) {
      y <- .subset2(.subset2(unit, gg), "eps");
      y <- log2(y);
      mad <- 1.4826 * median(abs(y));
      if (mad==0) {
        return(matrix(data=1, nrow=nrow(y), ncol=ncol(y)));
      }
      if (mergeGroups) {
        return(matrix(MASS_psi.huber(y/madMerged), ncol=ncol(y)));
      } else {
        return(matrix(MASS_psi.huber(y/mad), ncol=ncol(y)));
      }
    })
    res;
  } # resFcn()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  rs <- calculateResidualSet(this, verbose=verbose);
  ws <- getWeightsSet(this, verbose=verbose);
  nbrOfArrays <- length(rs);

  verbose && enter(verbose, "Calculating PLM weights");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get data and parameter objects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  if (is.null(ds)) {
    throw("No data set specified for PLM: ", getFullName(this));
  }

  cdf <- getCdf(ds);
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
  } else {
    nbrOfUnits <- length(units);
  }

  if (force) {
    unitsToDo <- units;
  } else {
    unitsToDo <- findUnitsTodo(ws, units=units);
  }

  ## Already done?
  if (length(unitsToDo) == 0) {
    verbose && cat(verbose, "All weights already calculated. Skipping.")
    verbose && exit(verbose)
    return(invisible(ws))
  }

  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

  unitsPerChunk <- ram * 100000/length(getDataSet(this));
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n",
                    nbrOfChunks, unitsPerChunk);
  head <- 1:unitsPerChunk;

  verbose && enter(verbose, "Extracting unit data");
  count <- 1;
  while (length(unitsToDo) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", count, nbrOfChunks));

    # Identify units to process in this chunk
    if (length(unitsToDo) < unitsPerChunk) {
      head <- 1:length(unitsToDo);
    }
    units <- unitsToDo[head];
    verbose && cat(verbose, "Number of units: ", length(units));

    residualsList <- readUnits(rs, units=units, verbose=less(verbose), stratifyBy="pm");

    verbose && enter(verbose, "Calculating weights");
    weightsList <- lapply(residualsList, FUN=resFcn, mergeGroups=this$mergeGroups);
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing weights");

    cdf <- getCellIndices(getCdf(ds), units=units, stratifyBy="pm", ...);

    for (ii in seq_along(ds)) {
      wf <- ws[[ii]];

      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", ii, getName(wf), length(ds)));

      data <- lapply(weightsList, function(unit) {
        lapply(unit, function(group) {
          nrow <- nrow(group);
          list(
            intensities=2^group[,ii],
            stdvs=rep(1, times=nrow),
            pixels=rep(1, times=nrow)
          );
        });
      });

      .updateCelUnits(getPathname(wf), cdf=cdf, data=data);

      verbose && exit(verbose);
    } # for (ii ...)

    verbose && exit(verbose);

    unitsToDo <- unitsToDo[-head];
    count <- count + 1;

    verbose && exit(verbose);
  } # while (...)

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  ## Generate checksum files
  wsZ <- getChecksumFileSet(ws)

  verbose && exit(verbose);

  invisible(ws);
}, protected=TRUE)


##########################################################################
# HISTORY:
# 2011-03-01 [HB]
# o Harmonized the verbose output.
# 2007-02-15
# o Based on ProbeLevelModel.calculateResiduals
#   and QualityAssessmentModel.getWeights
##########################################################################
