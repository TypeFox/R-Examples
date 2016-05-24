setMethodS3("getSubsetOfUnits", "AffymetrixCdfFile", function(this, units=NULL, unitTypes=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  useGenomeInformation <- FALSE;
  if (is.null(units)) {
  } else if (is.character(units)) {
    units <- Arguments$getCharacter(units);
    if (units %in% c("-X", "-Y", "-XY")) {
      useGenomeInformation <- TRUE;
    } else {
      throw("Unknown value on argument 'units': ", units);
    }
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }

  # Argument 'unitTypes':
  if (!is.null(unitTypes)) {
    unitTypes <- Arguments$getIndices(unitTypes, range=c(0,99));
    unitTypes <- sort(unitTypes);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Identifying subset of units");

  verbose && cat(verbose, "Argument 'units':");
  verbose && str(verbose, units);
  verbose && cat(verbose, "Argument 'unitTypes': ", paste(unitTypes, collapse=", "));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get genome information annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (useGenomeInformation) {
    verbose && enter(verbose, "Getting GenomeInformation file");
    # Get the genome information (throws an exception if missing)
    gi <- getGenomeInformation(this);
#    verbose && print(verbose, gi);
    giChecksum <- getChecksum(gi);
    verbose && exit(verbose);
  } else {
    giChecksum <- NULL;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify subset of units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying units to include and exclude");

  # Include all units and exclude none by default
  unitsIncl <- NULL;
  unitsExcl <- NULL;

  verbose && cat(verbose, "Argument 'units':");
  if (is.null(units)) {
    verbose && str(verbose, units);
  } else if (is.character(units)) {
    verbose && cat(verbose, units);
    # Select by chromosome(s)?
    if (units %in% c("-X", "-Y", "-XY")) {
      verbose && enter(verbose, "Selecting units by genomic location");

      # Identify chromosomes to be excluded
      parts <- gsub("-", "", units);
      parts <- strsplit(parts, split="", fixed=TRUE)[[1]];
      parts <- unique(parts);
      chromosomes <- c("X"=23, "Y"=24, "M"=25)[parts];
      if (anyNA(chromosomes)) {
        throw("Unknown chromosomes: ", parts[is.na(chromosomes)]);
      }
      chromosomes <- sort(chromosomes);

      verbose && cat(verbose, "Chromosomes to exclude:");
      verbose && str(verbose, chromosomes);

      unitsExcl <- getUnitsOnChromosomes(gi, chromosomes, .checkArgs=FALSE);
      # Not needed anymore
      chromosomes <- NULL;

      verbose && exit(verbose);
    } else {
      throw("Internal error. This statement should never be reached.");
    }
  } else {
    verbose && str(verbose, units);
    unitsIncl <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }

  if (is.null(unitsIncl))
    unitsIncl <- seq_len(nbrOfUnits(this));

  verbose && cat(verbose, "Units to include:");
  verbose && str(verbose, unitsIncl);

  verbose && cat(verbose, "Units to exclude:");
  verbose && str(verbose, unitsExcl);

  units <- setdiff(unitsIncl, unitsExcl);
  # Not needed anymore
  unitsIncl <- unitsExcl <- NULL;

  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify units by type
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(unitTypes)) {
    verbose && enter(verbose, "Filtering by unit type");
    verbose && cat(verbose, "Chip type: ", getChipType(this));
    verbose && cat(verbose, "Unit types: ", paste(unitTypes, collapse=", "));
    verbose && cat(verbose, "Units interrogated:");
    verbose && str(verbose, units);
    ut <- getUnitTypes(this, units=units, verbose=less(verbose, 5));
    verbose && cat(verbose, "Unit types:");
    verbose && print(verbose, table(ut));
    verbose && str(verbose, ut);

    verbose && cat(verbose, "Keeping units of interest:");
    keep <- which(ut %in% unitTypes);
    # Not needed anymore
    ut <- NULL;
    units <- units[keep];
    # Not needed anymore
    keep <- NULL;
    verbose && str(verbose, units);

    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  units;
}, protected=TRUE)  # getSubsetOfUnits()


############################################################################
# HISTORY:
# 2008-07-26
# o No need to cache results; fast enough.
# o Added support to filter by unit types.
# o Created from getSubsetOfCellIndices().
############################################################################
