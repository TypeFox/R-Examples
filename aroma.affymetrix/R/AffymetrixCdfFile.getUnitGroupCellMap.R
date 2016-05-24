setMethodS3("getUnitGroupCellMap", "AffymetrixCdfFile", function(this, units=NULL, retNames=FALSE, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  flattenCellIndices <- function(cells, ..., retNames=FALSE, verbose=FALSE) {
    # Returning indices or names?
    if (!retNames) {
      verbose && enter(verbose, "Renaming group names to group indices");
      cells <- lapply(cells, FUN=function(unit) {
        groups <- .subset2(unit, 1);
        names(groups) <- seq_len(length(groups));
        list(groups=groups);
      });

      verbose && print(verbose, gc);
      verbose && exit(verbose);
    }

    # Flatten cell data
    verbose && enter(verbose, "Flattening cell data");
    cells <- unlist(cells, use.names=TRUE);
    names <- names(cells);
    names(cells) <- NULL;
    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);

    verbose && enter(verbose, "Extract unit and group names");
    # Do some tricks to clean up the names
    names <- gsub("([.]groups|indices*)", "", names);
    pattern <- "^(.*)[.](.*)[.](.*)$";

    units <- gsub(pattern, "\\1", names);
    groups <- gsub(pattern, "\\2", names);
    # Not needed anymore
    pattern <- names <- NULL; # Not needed anymore
    verbose && exit(verbose);

    if (!retNames) {
      verbose && enter(verbose, "Converting to indices");
      units <- as.integer(units);
      groups <- as.integer(groups);
    }
    verbose && exit(verbose);

    # Return data
    map <- data.frame(unit=units, group=groups, cell=cells);
    class(map) <- c("UnitGroupCellMap", class(map));

    map;
  } # flattenCellIndices()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'retNames':
  retNames <- Arguments$getLogical(retNames);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting (unit, group, cell) map");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this);

  # Look for results in file cache
  verbose && enter(verbose, "Checking cache");
  key <- list(method="getUnitGroupCellMap", class=class(this)[1],
                   chipType=chipType, units=units, retNames=retNames, ...);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getUnitGroupCellMap", chipType=chipType, units=units, retNames=retNames, ...);
  }
  dirs <- c("aroma.affymetrix", chipType);
  map <- NULL;
  if (!force)
     map <- loadCache(key=key, dirs=dirs);
  if (is.null(map)) {
    verbose && exit(verbose, suffix="...miss");
  } else {
    verbose && printf(verbose, "RAM: %.2fMB\n", object.size(map)/1024^2);
    verbose && exit(verbose, suffix="...hit");
  }

  # Not in cache?
  if (is.null(map)) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Retrieve matrix of (unit, group, cell) indices
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cells <- getCellIndices(this, units=units, ..., verbose=less(verbose));
    nbrOfUnits <- length(cells);
    verbose && printf(verbose, "Read %d units\n", nbrOfUnits);

    if (!retNames) {
      # Convert unit names to unit indices
      if (is.null(units))
        units <- seq_len(nbrOfUnits);
      names(cells) <- units;

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && enter(verbose, "Flattening cell indices to create cell map");
    map <- flattenCellIndices(cells, retNames=retNames, verbose=less(verbose));
    verbose && str(verbose, map);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Save to cache
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Save only results > 50kB
    if (object.size(map) > 50e3) {
      saveCache(map, key=key, dirs=dirs);
      verbose && cat(verbose, "Saved to file cache");
    }
  } # if (is.null(map))

  # Extract subset of units
##  if (!is.null(units)) {
##    mapUnits <- map[,"unit"];
##    if (retNames) {
##      # Convert unit names in the map to unit indices
##      allUnitNames <- getUnitNames(this);
##      mapUnits <- match(mapUnits, allUnitNames);
##    }
##    keep <- (mapUnits %in% units);
##    map <- map[keep,,drop=FALSE];
##  }

  verbose && printf(verbose, "RAM: %.2fMB\n", object.size(map)/1024^2);
  verbose && exit(verbose);

  map;
}, protected=TRUE)  # getUnitGroupCellMap()



setMethodS3("getUnitGroupCellChromosomePositionMap", "AffymetrixCdfFile", function(this, units=NULL, chromosomes=NULL, orderByPosition=TRUE, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  ugcMap <- NULL;
  if (is.null(units)) {
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- ugcMap[,"unit"];
  }
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }
  units0 <- units;

  # Get the genome position information
  gi <- getGenomeInformation(this);

  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
    allChromosomes <- getChromosomes(gi);
    unknown <- chromosomes[!(chromosomes %in% allChromosomes)];
    if (length(unknown) > 0) {
      throw("Argument 'chromosomes' contains unknown values: ",
                                 paste(unknown, collapse=", "));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting (unit, group, cell, chromosome, position) map");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Look for results in file cache
  verbose && enter(verbose, "Checking cache");
  chipType <- getChipType(this);
  key <- list(method="getUnitGroupCellChromosomePositionMap",
              class=class(this)[1],
              chipType=chipType, units=units, ugcMap=ugcMap,
              chromosomes=chromosomes, orderByPosition=orderByPosition);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getUnitGroupCellChromosomePositionMap", chipType=chipType, units=units, ugcMap=ugcMap, chromosomes=chromosomes, orderByPosition=orderByPosition);
  }
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    map <- loadCache(key=key, dirs=dirs);
    if (!is.null(map)) {
      verbose && cat(verbose, "Found cached results");
      verbose && exit(verbose);
      return(map);
    }
  }


  # Select by chromosome(s)?
  if (!is.null(chromosomes)) {
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);
    verbose && cat(verbose, "Subset by chromosomes:");
    verbose && str(verbose, chromosomes);
    units <- getUnitsOnChromosomes(gi, chromosomes);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);
    if (!is.null(units0)) {
      units <- intersect(units, units0);
    }
  }
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);


  # Get the (unit, group, cell) map?
  if (isUnitGroupCellMap(ugcMap)) {
    ugcMap <- getUnitGroupCellMap(this, units=units, force=force, verbose=less(verbose, 10));
    verbose && cat(verbose, "(unit, group, cell) map:");
    verbose && str(verbose, ugcMap);
  }

  # Get the (chromosome, position) map
  cpMap <- getData(gi, units=ugcMap[,"unit"], force=force,
                                              verbose=less(verbose, 10));
  verbose && cat(verbose, "(chromosome, position) map:");
  verbose && str(verbose, cpMap);

  # Sanity check
  stopifnot(nrow(ugcMap) == nrow(cpMap));

  # Merge the two maps
  map <- cbind(ugcMap, cpMap);
  # Not needed anymore
  ugcMap <- cpMap <- NULL;

  if (orderByPosition) {
    o <- with(map, order(chromosome, physicalPosition));
    map <- map[o,,drop=FALSE];
    # Not needed anymore
    o <- NULL;
    verbose && cat(verbose, "Reordered by genomic position");
  }
  rownames(map) <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save only results > 50kB
  if (object.size(map) > 50e3) {
    saveCache(map, key=key, dirs=dirs);
    verbose && cat(verbose, "Saved to file cache");
  }

  verbose && exit(verbose);

  map;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-04-21
# o Now getUnitGroupCellMap() always returns a data.frame.
# o Added class attribute to getUnitGroupCellMap().
# 2008-03-11
# o Added getUnitGroupCellChromosomePositionMap();
# 2007-03-08
# o Note, getUnitGroupCellMap() is just a test function.
# o Created.
############################################################################
