setMethodS3("writeCdfByExcludingCells", "AffymetrixCdfFile", function(this, tags=c("filtered"), cellsToExclude, dropEmptyGroups=TRUE, dropEmptyUnits=dropEmptyGroups, ..., overwrite=FALSE, ram=1, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   maxNbrOfCells <- nbrOfCells(this);

  # Argument 'tags':
  tags <- Arguments$getCharacters(tags);
  tags <- tags[nzchar(tags)];

  # Argument 'cellsToExclude':
  cellsToExclude <- Arguments$getIndices(cellsToExclude, max=maxNbrOfCells);
  cellsToExclude <- unique(cellsToExclude);
  cellsToExclude <- sort(cellsToExclude);

  # Argument 'dropEmptyGroups':
  dropEmptyGroups <- Arguments$getLogical(dropEmptyGroups);

  # Argument 'dropEmptyUnits':
  dropEmptyUnits <- Arguments$getLogical(dropEmptyUnits);

  # Argument 'ram':
  ram <- Arguments$getDouble(ram, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Creating new CDF by dropping certain cells from existing CDF");
  chipType <- getChipType(this, fullname=TRUE);
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Added tags: ", paste(tags, collapse=","));

  pathname <- getPathname(this);
  verbose && cat(verbose, "Source CDF: ", getFilename(this));

  path <- getPath(this);
  fullnameF <- paste(c(chipType, tags), collapse=",");
  filenameF <- sprintf("%s.cdf", fullnameF);
  pathnameF <- Arguments$getWritablePathname(filenameF, path=path, mustNotExist=!overwrite);
  verbose && cat(verbose, "Destination CDF: ", filenameF);
  verbose && cat(verbose, "Path: ", path);

  # Get CDF header (to be reused when writing the new CDF)
  cdfHeader <- .readCdfHeader(pathname);

  # Get CDF QC units (to be reused when writing the new CDF)
  cdfQcUnits <- .readCdfQc(pathname);

  verbose2 <- as.logical(verbose);

  nbrOfUnits <- nbrOfUnits(this);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  verbose && cat(verbose, "Number of QC units (not filtered): ", length(cdfQcUnits));
  cells <- getCellIndices(this, unlist=TRUE, useNames=FALSE);
  nbrOfCells <- length(cells);
  verbose && cat(verbose, "Number of cells in input CDF: ", nbrOfCells);
  verbose && cat(verbose, "Requestes cells to be excluded:");
  verbose && str(verbose, cellsToExclude);

  cellsToExclude <- intersect(cellsToExclude, cells);
  verbose && cat(verbose, "Requestes cells to be excluded in input CDF:");
  verbose && str(verbose, cellsToExclude);

  verbose && printf(verbose, "Number of cells to exclude: %d (%.1f%%) out of %d\n", length(cellsToExclude), 100*length(cellsToExclude)/nbrOfCells, nbrOfCells);

  nbrOfUnitsPerChunk <- ram*100e3;
  units <- seq_len(nbrOfUnits);
  nbrOfChunks <- ceiling(nbrOfUnits / nbrOfUnitsPerChunk);
  nbrOfUnits <- nbrOfUnits(this);
  verbose && cat(verbose, "Number of units per chunk: ", nbrOfUnitsPerChunk);
  head <- 1:nbrOfUnitsPerChunk;

  # Allocating output CDF list structure
  cdfListF <- vector("list", length=nbrOfUnits);
  names(cdfListF) <- getUnitNames(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading CDF, dropping cells, ...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfUnitsLeft <- length(units);
  chunk <- 1L;
  while (nbrOfUnitsLeft > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", chunk, nbrOfChunks));

    # Special case for last chunk
    if (nbrOfUnitsLeft < nbrOfUnitsPerChunk) {
      head <- 1:nbrOfUnitsLeft;
    }

    unitsChunk <- units[head];
    units <- units[-head];
    nbrOfUnitsLeft <- length(units);

    verbose && cat(verbose, "Units:");
    verbose && str(verbose, unitsChunk);

    # Sanity check
    stopifnot(length(unitsChunk) > 0);

    verbose && enter(verbose, "Reading subset of units from source CDF");
    cdfListChunk <- .readCdf(pathname, readIndices=TRUE, units=unitsChunk);
    verbose && exit(verbose);

    verbose && enter(verbose, "Dropping cells");
    cdfListChunkF <- dropCellsFromCdfList(cdfListChunk, cellsToExclude=cellsToExclude, maxNbrOfCells=maxNbrOfCells, verbose=verbose2);
    # Not needed anymore
    cdfListChunk <- NULL; # Not needed anymore
    verbose && exit(verbose);

    # Drop empty groups?
    if (dropEmptyGroups) {
      verbose && enter(verbose, "Dropping empty unit groups");
      nbrOfGroups <- 0L;
      nbrDropped <- 0L;
      for (uu in seq_along(cdfListChunkF)) {
        unit <- cdfListChunkF[[uu]];
        groups <- unit$groups;
        nbrOfGroups <- nbrOfGroups + length(groups);
        ns <- lapply(groups, FUN=function(group) length(group$indexpos));
        groups <- groups[(ns > 0)];
        unit$groups <- groups;
        cdfListChunkF[[uu]] <- unit;
        nbrDropped <- nbrDropped + sum(ns == 0);
      } # for (uu ...)
      verbose && printf(verbose, "Number of empty unit groups dropped: %d (%.1f%%) of %d\n", nbrDropped, 100*nbrDropped/nbrOfGroups, nbrOfGroups);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Appending to result CDF structure");
    cdfListF[unitsChunk] <- cdfListChunkF;
    # Not needed anymore
    cdfListChunkF <- NULL; # Not needed anymore
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    chunk <- chunk + 1L;

    verbose && exit(verbose);
  } # while (nbrOfUnitsLeft > 0)

  # Sanity check
  stopifnot(length(cdfListF) == nbrOfUnits);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop empty units?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (dropEmptyUnits) {
    verbose && enter(verbose, "Dropping empty units");

    ns <- lapply(cdfListF, FUN=function(unit) {
      length(unit$groups);
    });
    cdfListF <- cdfListF[(ns > 0)];
    nbrOfUnitsF <- length(cdfListF);

    nbrDropped <- nbrOfUnits-nbrOfUnitsF;
    verbose && printf(verbose, "Number of empty units dropped: %d (%.1f%%) of %d\n", nbrDropped, 100*nbrDropped/nbrOfUnits, nbrOfUnits);
    verbose && cat(verbose, "Number of remaining units: ", nbrOfUnitsF);

    # Update CDF header
    cdfHeader$nunits <- nbrOfUnitsF;
    cdfHeader$probesets <- nbrOfUnitsF;

    verbose && exit(verbose);
  } # if (dropEmptyUnits)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Writing CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing filtered CDF");
  if (overwrite && isFile(pathnameF)) {
     file.remove(pathnameF);
  }

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathnameF, verbose=verbose);

  .writeCdf(pathnameT, cdfheader=cdfHeader, cdf=cdfListF, cdfqc=cdfQcUnits, verbose=verbose2);

  # Not needed anymore
  cdfHeader <- cdfListF <- NULL; # Not needed anymore

  # Validate new CDF
  cdfT <- newInstance(this, pathnameT);
  cellsF <- getCellIndices(cdfT, unlist=TRUE, useNames=FALSE);

  # Sanity check
  stopifnot(!any(is.element(cellsToExclude, cellsF)));

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=verbose);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return new CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdfF <- newInstance(this, pathnameF);
  verbose && print(verbose, cdfF);

  ## Create checksum file
  cdfFZ <- getChecksumFile(cdfF)

  nbrOfCellsF <- length(cellsF);
  verbose && cat(verbose, "Number of cells in output CDF: %d", nbrOfCellsF);
  verbose && printf(verbose, "Number of cells excluded: %d (%.1f%%) of %d\n", nbrOfCells-nbrOfCellsF, 100*(nbrOfCells-nbrOfCellsF)/nbrOfCells, nbrOfCells);

  verbose && exit(verbose);

  cdfF;
}, protected=TRUE) # writeCdfByExcludingCells()


############################################################################
# HISTORY:
# 2011-09-17
# o Now handling empty groups and empty units as well.
# o Added writeCdfByExcludingCells() for AffymetrixCdfFile.
# o Created.
############################################################################
