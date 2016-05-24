setMethodS3("getUnitGroupNamesFromUgcMap", "AffymetrixCdfFile", function(this, ugcMap, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting unit and group names from CDF");
  allUnits <- ugcMap[,"unit"];

  # Get unit names
  res <- data.frame(
    unitName = getUnitNames(this, units=allUnits),
    groupName = character(length(allUnits)),
    stringsAsFactors = FALSE
  );
  verbose && str(verbose, res);

  verbose && enter(verbose, "Reading all group names for units of interest");
  uniqueUnits <- unique(allUnits);
  groupNames <- .readCdfGroupNames(getPathname(this), units=uniqueUnits);
  verbose && cat(verbose, "First unit:");
  verbose && str(verbose, groupNames[1]);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating the number of groups per unit");
  unitSizes <- lapply(groupNames, FUN=length);
  unitSizes <- unlist(unitSizes, use.names=FALSE);
  uniqueUnitSizes <- sort(unique(unitSizes));
  verbose && exit(verbose);

  verbose && enter(verbose, "Building (unit name, group name) map");

  for (uu in seq_along(uniqueUnitSizes)) {
    unitSize <- uniqueUnitSizes[uu];
    verbose && enter(verbose, "Processing set #%d of %d containing units with a maximum of %d group(s)", uu, length(uniqueUnitSizes), unitSize);

    # Extract the group names for unit with 'unitSize' groups as a matrix
    idxs <- which(unitSizes == unitSize);
    units <- uniqueUnits[idxs];
    verbose && cat(verbose, "Number of units with ", unitSize,
                                            " group(s): ", length(units));

    if (length(units) > 0) {  # Isn't this always the case? /HB 2008-04-28
      names <- groupNames[idxs];
      names <- unlist(names, use.names=FALSE);
      names <- matrix(names, nrow=unitSize);
      verbose && cat(verbose, "Identfied group names:");
      verbose && str(verbose, names);

      # Find the subset of the UGC map that contains these units
      rrU <- (ugcMap[,"unit"] %in% units);
#      verbose && str(verbose, rrU);

      # For each possible group index...
      for (group in seq_len(unitSize)) {
        verbose && enter(verbose, "Group ", group);
        # Identify row in the UGC map containing those units and the group
        rrG <- (ugcMap[,"group"] == group);
#        verbose && str(verbose, rrG);
        rrUG <- which(rrU & rrG);
#        verbose && str(verbose, rrUG);
        if (length(rrUG) > 0) {
          unitsUG <- ugcMap[rrUG, "unit"];
          rr <- match(unitsUG, units);
          res[rrUG,"groupName"] <- names[group,rr];
        }
        verbose && exit(verbose);
      }
    }

    verbose && exit(verbose);
  } # for (uu ...)

  verbose && exit(verbose);

  if (nrow(res) != nrow(ugcMap)) {
    throw("Internal error: Number of extract unit and group names does not match the number of rows in the UGC map: ", nrow(res), " != ", nrow(ugcMap));
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE) # getUnitGroupNamesFromUgcMap()


############################################################################
# HISTORY:
# 2011-11-18
# o Now the verbose progress output of getUnitGroupNamesFromUgcMap()
#   is more informative on what subset in order is currently processed.
# 2008-04-28
# o SPEEDUP: Changed the algorithm for getUnitGroupNamesFromUgcMap(). It
#   was painfully slow for large UGC maps.  Took ~10-14 days(!) for the
#   GenomeWideSNP_6 chip. Now 50s. The relative overhead for loading/saving
#   to cache is now so large that it no longer pays off.
# 2008-02-27
# o Now getUnitGroupNamesFromUgcMap() caches results to file.
# 2008-02-05
# o Created.
############################################################################
