setMethodS3("getUnitGroupCellMatrixMap", "ChipEffectFile", function(this, units=NULL, groups=NULL, ..., verbose=FALSE) {
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
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'groups':
  if (!is.null(groups)) {
    groups <- Arguments$getIndices(groups, max=999);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }


  # Subset by groups?
  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Build integer UxG matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  allUnits <- unique(ugcMap[,"unit"]);
  nbrOfGroups <- length(groups);
  naValue <- as.integer(NA);
  map <- matrix(naValue, nrow=nbrOfUnits, ncol=nbrOfGroups);
  
  for (gg in seq_len(nbrOfGroups)) {
    group <- groups[gg];
    verbose && enter(verbose, sprintf("Group %d (%d) of %d", 
                                                  gg, group, nbrOfGroups));

    idxs <- which(ugcMap$group == group);
    units <- ugcMap[idxs, "unit"];
    cells <- ugcMap[idxs, "cell"];
    rr <- match(units, allUnits);
    map[rr,gg] <- cells;

    # Not needed anymore
    idxs <- rr <- units <- cells <- NULL;
    verbose && exit(verbose);
  }

  class(map) <- "UnitGroupCellMatrixMap";

  verbose && cat(verbose, "Unit-by-group cell matrix map:");
  verbose && str(verbose, map);

  map;
}, protected=TRUE)  # getUnitGroupCellMatrixMap()



setMethodS3("getUnitGroupCellArrayMap", "ChipEffectFile", function(this, units=NULL, groups=NULL, groupCells=NULL, ..., verbose=FALSE) {
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
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'groups':
  if (!is.null(groups)) {
    groups <- Arguments$getIndices(groups, max=999);
  }

  # Argument 'groupCells':
  if (!is.null(groupCells)) {
    groupCells <- Arguments$getIndices(groupCells);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }


  # Subset by groups?
  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  # Subset by group cells (cells indexed within group)?
  if (!is.null(groupCells)) {
    idxs <- which(ugcMap$cell %in% groupCells);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groupCells <- sort(unique(ugcMap$cell));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Build integer UxG matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  allUnits <- unique(ugcMap[,"unit"]);
  nbrOfGroups <- length(groups);
  naValue <- as.integer(NA);
  map <- matrix(naValue, nrow=nbrOfUnits, ncol=nbrOfGroups);
  
  for (gg in seq_len(nbrOfGroups)) {
    group <- groups[gg];
    verbose && enter(verbose, sprintf("Group %d (%d) of %d", 
                                                  gg, group, nbrOfGroups));

    idxs <- which(ugcMap$group == group);
    units <- ugcMap[idxs, "unit"];
    cells <- ugcMap[idxs, "cell"];
    rr <- match(units, allUnits);
    map[rr,gg] <- cells;

    # Not needed anymore
    idxs <- rr <- units <- cells <- NULL;
    verbose && exit(verbose);
  }

  class(map) <- "UnitGroupCellMatrixMap";

  verbose && cat(verbose, "Unit-by-group cell matrix map:");
  verbose && str(verbose, map);

  map;
}, protected=TRUE)  # getUnitGroupCellArrayMap()



############################################################################
# HISTORY:
# 2008-07-13
# o Added argument 'drop=FALSE' to extractTheta().
# 2008-06-09
# o Added getUnitGroupCellMatrixMap() to ChipEffectFile.  The extractTheta()
#   methods is now using this method.
# 2008-05-10
# o Updated to take an UGC map via argument 'units'.
# 2008-05-09
# o Created.
############################################################################
