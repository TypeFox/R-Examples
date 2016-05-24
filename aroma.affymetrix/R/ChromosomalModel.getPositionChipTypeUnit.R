setMethodS3("getPositionChipTypeUnit", "ChromosomalModel", function(this, chromosome, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'chromosome':
  chromosome <- Arguments$getIndex(chromosome);
  knownChromosomes <- getChromosomes(this);
  if (!chromosome %in% knownChromosomes) {
    throw("Argument 'chromosome' contains an unknown chromosome: ", chromosome);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (position, chipType, unit) map for this chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting (position, chipType, unit) map");

  # Get the UnitNameFile:s
  unfList <- getListOfUnitNamesFiles(this, verbose=less(verbose, 10));

  # Get the genome information files
  ugpList <- lapply(unfList, FUN=getAromaUgpFile, verbose=less(verbose, 10));
  verbose && print(verbose, ugpList);

  # Get the genome information files
  ugpList <- lapply(unfList, FUN=getAromaUgpFile, verbose=less(verbose, 10));
  verbose && print(verbose, ugpList);

  # Get the units on the chromosome of interest
  unitsList <- lapply(ugpList, FUN=function(ugp) {
    getUnitsOnChromosome(ugp, chromosome=chromosome, ...);
  });
  verbose && str(verbose, unitsList);

  # Gets (position, chipType) for these units
  posList <- vector("list", length(unitsList));
  names(posList) <- names(unitsList);
  chipTypeList <- vector("list", length(unitsList));
  names(chipTypeList) <- names(unitsList);
  for (kk in seq_along(posList)) {
    ugp <- ugpList[[kk]];
    units <- unitsList[[kk]];
    pos <- getPositions(ugp, units=units);

    # Keep only units with a position
    keep <- which(is.finite(pos));
    nbrOfUnitsBefore <- length(pos);
    nbrOfUnits <- length(keep);
    nbrOfUnitsExcl <- nbrOfUnitsBefore - nbrOfUnits;
    if (nbrOfUnitsExcl > 0) {
      pos <- pos[keep];
      units <- units[keep];
      verbose && cat(verbose, "Excluded ", nbrOfUnitsExcl, " (out of", nbrOfUnitsBefore, ") units because there is no position information available for those.");
    }
    unitsList[[kk]] <- units;
    posList[[kk]] <- pos;
    chipTypeList[[kk]] <- rep(kk, length(units));
    # Not needed anymore
    ugp <- units <- keep <- NULL;
  }
  # Not needed anymore
  ugpList <- NULL;

  verbose && str(verbose, unitsList);
  verbose && str(verbose, posList);
  verbose && str(verbose, chipTypeList);

  # Unlist and order (units, position, chipType) by position
  pos <- unlist(posList, use.names=FALSE);
  # Not needed anymore
  posList <- NULL;
  o <- order(pos);
  pos <- pos[o];

  chipType <- unlist(chipTypeList, use.names=FALSE);
  # Not needed anymore
  chipTypeList <- NULL;
  chipType <- chipType[o];

  # Convert chipType into a factor
  chipTypes <- sapply(unfList, FUN=getChipType);
  attr(chipType, "levels") <- chipTypes;
  class(chipType) <- "factor";
  # Not needed anymore
  unfList <- NULL;

  units <- unlist(unitsList, use.names=FALSE);
  # Not needed anymore
  unitsList <- NULL;
  units <- units[o];
  # Not needed anymore
  o <- NULL;

  pcu <- data.frame(position=pos, chipType=chipType, unit=units);
  # Not needed anymore
  units <- pos <- chipType <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && cat(verbose, "(position, chipType, unit) map:");
  verbose && str(verbose, pcu);

  verbose && exit(verbose);

  pcu;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2009-01-26
# o Updated getPositionChipTypeUnit() of ChromosomalModel to utilize
#   the UnitNamesFile Interface instead of assuming an AffymetrixCdfFile.
#   This requires aroma.core v1.0.1.
# 2007-09-25
# o Moved getPositionChipTypeUnit() to ChromosomalModel.
# 2007-09-20
# o Added getPositionChipTypeUnit().
############################################################################
