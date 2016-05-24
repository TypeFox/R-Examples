setMethodS3("getSubsetOfCellIndices", "AffymetrixCdfFile", function(this, units=NULL, stratifyBy=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
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

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Identifying subset of cell indices");

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
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="getSubsetOfCellIndices", class=class(this)[1],
              chipType=getChipType(this), giChecksum=giChecksum,
              units=units, stratifyBy=stratifyBy, ...);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getSubsetOfCellIndices", chipType=getChipType(this), giChecksum=giChecksum, units=units, stratifyBy=stratifyBy, ...);
  }
  dirs <- c("aroma.affymetrix", getChipType(this));
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Found cached results.");
      verbose && exit(verbose);
      return(res);
    }
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

  verbose && cat(verbose, "Units to include:");
  verbose && str(verbose, unitsIncl);

  verbose && cat(verbose, "Units to exclude:");
  verbose && str(verbose, unitsExcl);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying cell indices to include and exclude
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying cells to include and exclude");

  verbose && cat(verbose, "Stratify by cell type(s):");
  verbose && str(verbose, stratifyBy);

  if (is.null(unitsIncl) & is.null(unitsExcl)) {
    # All cells by default
    cellsIncl <- NULL;
    cellsExcl <- NULL;
  } else if (is.null(unitsIncl) & !is.null(unitsExcl)) {
    # All cells but...
    cellsIncl <- NULL;
    verbose && enter(verbose, "Reading cell indices to exclude");
    cellsExcl <- getCellIndices(this, units=unitsExcl,
                   stratifyBy=stratifyBy, useNames=FALSE, unlist=TRUE);
    verbose && exit(verbose);
  } else if (!is.null(unitsIncl) & is.null(unitsExcl)) {
    verbose && enter(verbose, "Reading cell indices to include");
    cellsIncl <- getCellIndices(this, units=unitsIncl,
                   stratifyBy=stratifyBy, useNames=FALSE, unlist=TRUE);
    verbose && exit(verbose);
    cellsExcl <- NULL;
  } else if (!is.null(unitsIncl) & !is.null(unitsExcl)) {
    unitsIncl <- setdiff(unitsIncl, unitsExcl);
    verbose && enter(verbose, "Reading cell indices to include");
    cellsIncl <- getCellIndices(this, units=unitsIncl,
                   stratifyBy=stratifyBy, useNames=FALSE, unlist=TRUE);
    verbose && exit(verbose);
    cellsExcl <- NULL;
  }

  # Not needed anymore
  # Not needed anymore
  unitsIncl <- unitsExcl <- NULL;


  if (is.null(cellsIncl)) {
    # All types cells of cells?
    if (is.null(stratifyBy)) {
      cellsIncl <- seq_len(nbrOfCells(this));
    } else {
      verbose && enter(verbose, "Reading cell indices to include");
      cellsIncl <- getCellIndices(this, stratifyBy=stratifyBy,
                                            useNames=FALSE, unlist=TRUE);
      verbose && exit(verbose);
    }
  }
  verbose && cat(verbose, "Cells to include:");
  verbose && str(verbose, cellsIncl);

  verbose && cat(verbose, "Cells to exclude:");
  verbose && str(verbose, cellsExcl);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Including and excluding cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cells <- setdiff(cellsIncl, cellsExcl);
  # Not needed anymore
  cellsIncl <- cellsExcl <- NULL;

  verbose && cat(verbose, "Final set of cell indices:");
  verbose && str(verbose, cells);


  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    saveCache(cells, key=key, dirs=dirs);
  }

  cells;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-12-03
# o getSubsetOfCellIndices() stored in memoization in the root cache path.
# 2008-07-16
# o Created.
############################################################################
