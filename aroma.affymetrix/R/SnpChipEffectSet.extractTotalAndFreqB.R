setMethodS3("extractTotalAndFreqB", "CnChipEffectSet", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    ugcMap <- NULL;
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
    nbrOfUnits <- length(units);
    # Not needed anymore
    units <- NULL;
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting (total, freqB)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify possible groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (this$combineAlleles && this$mergeStrands) {
    # theta == (theta)
    groups <- 1;
  } else if (this$combineAlleles && !this$mergeStrands) {
    # theta == (theta+, theta-)
    groups <- c(1,3);
  } else if (!this$combineAlleles && this$mergeStrands) {
    # theta == (thetaA, thetaB)
    groups <- c(1,2);
  } else if (!this$combineAlleles && !this$mergeStrands) {
    # theta == (thetaA+, thetaB+, thetaA-, thetaB-)
    groups <- c(1,2,3,4);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(this);
  dim <- c(nbrOfUnits, 2, nbrOfArrays);
  dimnames <- list(NULL, c("total", "freqB"), getNames(this));
  data <- array(NA, dim=dim, dimnames=dimnames);

  for (kk in seq_len(nbrOfArrays)) {
    ce <- this[[kk]];
    dataKK <- extractTotalAndFreqB(ce, units=ugcMap, ...,
                                                verbose=less(verbose, 5));
    verbose && str(verbose, dataKK);
    data[,,kk] <- dataKK;
    # Not needed anymore
    dataKK <- NULL;
  }
  # Not needed anymore
  ugcMap <- NULL;

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }

  verbose && cat(verbose, "Results:");
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
})





setMethodS3("extractTotalAndFreqB", "SnpChipEffectSet", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    ugcMap <- NULL;
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
    nbrOfUnits <- length(units);
    # Not needed anymore
    units <- NULL;
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Extracting (total, freqB)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify possible groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (this$mergeStrands) {
    # theta == (thetaA, thetaB)
    groups <- c(1,2);
  } else {
    # theta == (thetaA+, thetaB+, thetaA-, thetaB-)
    groups <- c(1,2,3,4);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(this);
  dim <- c(nbrOfUnits, 2, nbrOfArrays);
  dimnames <- list(NULL, c("total", "freqB"), getNames(this));
  data <- array(NA, dim=dim, dimnames=dimnames);

  for (kk in seq_len(nbrOfArrays)) {
    ce <- this[[kk]];
    dataKK <- extractTotalAndFreqB(ce, units=ugcMap, ...,
                                                verbose=less(verbose, 5));
    verbose && str(verbose, dataKK);
    data[,,kk] <- dataKK;
    # Not needed anymore
    dataKK <- NULL;
  }
  # Not needed anymore
  ugcMap <- NULL;

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }

  verbose && cat(verbose, "Results:");
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
})



############################################################################
# HISTORY:
# 2008-12-09
# BUG FIX: There was no extractTotalAndFreqB() for CnChipEffectSet, but
# only for CnChipEffecSet (misspelled). Thanks Pierre Neuvial for spotting
# this.
# 2008-07-16
# o Added argument 'drop=FALSE' to all extractTotalAndFreqB().
# 2008-05-10
# o Created.
############################################################################
