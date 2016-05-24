###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod createUniqueCdf
# @aliasmethod createUnique
#
# @title "Creates a unique-cell version of the CDF"
#
# \description{
#  In some cases, single probes map to multiple genome locations.  In cases
#  where you may later want to store a probe estimate (e.g. a probe effect
#  or a residual), you will not be able to store more than one per cell.
#  The unique CDF facilitates this by making the cell indices unique
#  (essentially copying the multimapping probes).
# }
#
# @synopsis
#
# \arguments{
#   \item{chipType}{The chip type of the new CDF.}
#   \item{tags}{Tags added to the chip type of the new CDF.}
#   \item{path}{The path where to store the new CDF file.}
#   \item{...}{Additional arguments passed to @see "affxparser::writeCdf".}
#   \item{ram}{A @double saying if more or less units should be converted
#      per chunk.  A smaller value uses less memory.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the unique CDF as an @see "AffymetrixCdfFile" object.
# }
#
# \details{
#   Note that the set of units in the "unique" CDF is identical to that
#   of the input CDF.  So are the unit groups in each unit.
#   Also the number of cells per unit group is preserved.
#   It is only the cell-index composition of each unit group that changes.
#   The cell indices in the unique CDF are guaranteed to occur only
#   once, whereas this may not be true for the input CDF.
# }
#
# @author "MR"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("createUniqueCdf", "AffymetrixCdfFile", function(this, chipType=getChipType(this), tags="unique", path=NULL, units=NULL, ..., ram=NULL, verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rearrangeCells <- function(units, offset=0, hasGroups=TRUE, ncols, ...) {
    rearrangeGroup <- function(group, idxs, ...) {
      y = (idxs-1) %/% ncols;
      x = (idxs-1) - ncols*y;
      group$y <- y;
      group$x <- x;
      group$indices <- idxs;
      group;
    } # rearrangeGroup()

    nbrOfCells <- lapply(units, FUN=function(unit) .subset2(unit, "ncells"));
    nbrOfCells <- sum(unlist(nbrOfCells, use.names=FALSE));

    cells <- seq(from=offset+1, to=offset+nbrOfCells);

    verbose && printf(verbose, "Units: ");
    if (hasGroups) {
      for (kk in seq_along(units)) {
        if (verbose) {
          if (kk %% 1000 == 0) {
            printf(verbose, "%d, ", kk);
          } else if (kk %% 100 == 0) {
            cat(".");
          }
        }
        # groups <- units[[kk]]$groups;
        groups <- .subset2(.subset2(units, kk), "groups");
        for (ll in seq_along(groups)) {
          group <- .subset2(groups, ll);
          # Number of cells in this group
          nindices <- length(.subset2(group, "indices"));
          head <- 1:nindices;
          idxs <- .subset(cells, head);
          cells <- .subset(cells, (nindices+1):length(cells));
          groups[[ll]] <- rearrangeGroup(group, idxs);
        }
        units[[kk]]$groups <- groups;
      }
    } else {
      for (kk in seq_along(units)) {
        if (verbose) {
          if (kk %% 1000 == 0) {
            printf(verbose, "%d, ", kk);
          } else if (kk %% 100 == 0) {
            cat(".");
          }
        }
        group <- .subset2(units, kk);
        # Number of cells in this group
        nindices <- length(.subset2(group, "indices"));
        head <- 1:nindices;
        idxs <- .subset(cells, head);
        cells <- .subset(cells, (nindices+1):length(cells));
        group <- rearrangeGroup(group, idxs);
        units[[kk]] <- group;
      }
    }
    verbose && printf(verbose, "\n");

    units;
  } # rearrangeCells()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  verbose2 <- as.integer(verbose)-1;  # For 'affxparser' calls.
  verbose2 <- 2;

  # Argument 'units':
  if(is.null(units)) {
    units <- seq_len(nbrOfUnits(this));
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }



  # Already a unique CDF?
  if (isUniqueCdf(this)) {
    return(this);
  }

  verbose && enter(verbose, "Creating unique CDF");

  verbose && cat(verbose, "Chip type: ", getChipType(this));

  # Get the pathname of the source
  src <- getPathname(this);
  src <- Arguments$getReadablePathname(src);

  # Create the pathname of the destination CDF
  if (is.null(path)) {
    mainChipType <- gsub("[,].*", "", chipType);
    path <- filePath("annotationData", "chipTypes", mainChipType);
  }

  # Write to a temporary file first, and rename it when we know its complete
  name <- paste(c(chipType, tags), collapse=",");
  filename <- sprintf("%s.CDF", name);
  pathname <- Arguments$getWritablePathname(filename, path=path);
  pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=TRUE);

  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

  # Assure source and destination is not identical
  if (identical(src, pathnameT)) {
    throw("Cannot not create CDF file. Destination is same as source: ", src);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Figure out the number of (regular) unit cells in the output
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the number of groups per units
  verbose && enter(verbose, "Reading CDF group names");
  nbrOfGroupsPerUnit <- .readCdfGroupNames(src);
  verbose && exit(verbose);

  names(nbrOfGroupsPerUnit) <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Get the number of groups per unit
  nbrOfGroupsPerUnit <- sapply(nbrOfGroupsPerUnit, FUN=length);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # get a relatively low-memory version of the CDF indices
  # Q: Why is readCdf(..., stratifyBy=<vector>) not throwing?
  #    Error in match.arg(stratifyBy) : 'arg' must be of length 1
  #    /HB 2011-04-15
  cdfLite <- .readCdf(src, units=units,
              readXY=FALSE, readBases=FALSE,
              readIndexpos=FALSE, readAtoms=FALSE,
              readUnitType=FALSE, readUnitDirection=FALSE,
              readUnitNumber=FALSE, readUnitAtomNumbers=FALSE,
              readGroupAtomNumbers=FALSE, readGroupDirection=FALSE,
              readIndices=TRUE, readIsPm=FALSE,
              stratifyBy=c("nothing", "pmmm", "pm", "mm"), verbose=0);

  # Get the number of cells per unit
  nbrOfCellsPerUnit <- sapply(cdfLite, FUN=function(unit) {
    length(unlist(unit, use.names=FALSE));
  });
  verbose && cat(verbose, "Number of cells per unit:");
  verbose && summary(verbose, nbrOfCellsPerUnit);

  # Not needed anymore
  cdfLite <- NULL; # Not needed anymore

  # Total number of cells in CDF
  nbrOfCells <- sum(nbrOfCellsPerUnit);

  # Number of units in CDF
  nbrOfUnits <- length(nbrOfCellsPerUnit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # QC units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Keep all QC units
  verbose && enter(verbose, "Reading CDF QC units");
  destQcUnits <- .readCdfQc(src);
  verbose && exit(verbose);
  nbrOfQcUnits <- length(destQcUnits);
  nbrOfCellsPerQcUnit <- lapply(destQcUnits, FUN=.subset2, "ncells");
  nbrOfCellsPerQcUnit <- unlist(nbrOfCellsPerQcUnit, use.names=FALSE);
  nbrOfQcCells <- sum(nbrOfCellsPerQcUnit);
  verbose && printf(verbose,
                        "Number of QC cells: %d in %d QC units (%.1fMB)\n",
              nbrOfQcCells, nbrOfQcUnits, object.size(destQcUnits)/1024^2);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chip layout
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  totalNbrOfCells <- nbrOfCells + nbrOfQcCells;
  verbose && printf(verbose, "Total number of cells: %d\n", totalNbrOfCells);

  # Figure out a best fit square(!) layout of this number of cells
  side <- as.integer(floor(sqrt(totalNbrOfCells)));
  nrows <- ncols <- side;

  # If not big enough, increase the nbr of rows, then nbr of column
  if (nrows*ncols < totalNbrOfCells) {
    nrows <- as.integer(nrows + 1);
    if (nrows*ncols < totalNbrOfCells) {
      ncols <- as.integer(ncols + 1);
    }
  }
  verbose && printf(verbose, "Best array dimension: %dx%d (=%d cells, i.e. %d left-over cells)\n", nrows, ncols, nrows*ncols, nrows*ncols - totalNbrOfCells);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create CDF header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating CDF header with source CDF as template");

  verbose && enter(verbose, "Setting up header");

  # Get template header
  verbose && enter(verbose, "Reading CDF header");
  destHeader <- .readCdfHeader(src);
  verbose && exit(verbose);

  # Update CDF header dimension
  destHeader$nrows <- nrows;
  destHeader$ncols <- ncols;

  # Get unit names
  verbose && enter(verbose, "Reading CDF unit names");
  unitNames <- .readCdfUnitNames(src);
  verbose && exit(verbose);

  if (nbrOfUnits > 0) {
    # Number of bytes per unit:
    #  = 20 + (18+64)*nbrOfGroupsPerUnit + 14*nbrOfCellsPerUnit bytes
    unitLengths <- 20 + 82*nbrOfGroupsPerUnit + 14*nbrOfCellsPerUnit;
    avgUnitLength <- mean(unitLengths);
  } else {
    unitLengths <- NULL;
  }

  # Number of bytes per QC unit:
  if (nbrOfQcUnits > 0) {
    # Start positions for QC units
    # Number of bytes per QC unit:
    #  = 6 + 7*nbrOfCellsPerQcUnit bytes
    qcUnitLengths <- 6 + 7*nbrOfCellsPerQcUnit;
  } else {
    qcUnitLengths <- NULL;
  }
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing");
  verbose && cat(verbose, "destHeader:");
  verbose && str(verbose, destHeader);
  verbose && cat(verbose, "unitNames:");
  verbose && str(verbose, unitNames);
  verbose && cat(verbose, "qcUnitLengths:");
  verbose && str(verbose, qcUnitLengths);
  verbose && cat(verbose, "unitLengths:");
  verbose && str(verbose, unitLengths);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Open output connection
  con <- file(pathnameT, open = "wb");
  on.exit({
    if (!is.null(con))
      close(con);
    con <- NULL;
  });

  # Allocate/write CDF header
  .writeCdfHeader(con=con, destHeader, unitNames=unitNames,
                    qcUnitLengths=qcUnitLengths, unitLengths=unitLengths,
                                                        verbose=verbose2);
  # Not needed anymore
  # Not needed anymore
  destHeader <- unitNames <- qcUnitLengths <- unitLengths <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Restructure QC units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing QC units");
  verbose && enter(verbose, "Rearranging QC unit cell indices");
  destQcUnits <- rearrangeCells(destQcUnits, offset=nbrOfCells, ncols=ncols,
                                         hasGroups=FALSE, verbose=verbose);
  # Garbage collect
  gc <- gc();
  verbose && exit(verbose);

  # Write QC units
  .writeCdfQcUnits(con=con, destQcUnits, verbose=verbose2);
  # Not needed anymore
  destQcUnits <- NULL; # Not needed anymore

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create all new units in chunks
  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);
  verbose && printf(verbose, "Argument 'ram': %f\n", ram);
  verbose && printf(verbose, "Average unit length: %f bytes\n", avgUnitLength);

  # Scale (not used)
  if (nbrOfUnits > 0) {
    scale <- 200/avgUnitLength;
  } else {
    scale <- 1;
  }

  # Number of units per chunk
  nbrOfUnitsPerChunk <- as.integer(ram * scale * 20000);

  nbrOfChunks <- ceiling(nbrOfUnits / nbrOfUnitsPerChunk);
  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n",
                                         nbrOfChunks, nbrOfUnitsPerChunk);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading and extracting unit data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading, extracting, and writing units");

#  fields <- c("x", "y", "indices", "pbase", "tbase", "atom", "indexpos");
  fields <- c("pbase", "tbase", "atom", "indexpos");
  head <- 1:nbrOfUnitsPerChunk;
  count <- 1;
  idxOffset <- as.integer(0);
  unitsToDo <- 1:nbrOfUnits;
  while (length(unitsToDo) > 0) {
    if (length(unitsToDo) < nbrOfUnitsPerChunk) {
      head <- 1:length(unitsToDo);
    }
    units <- unitsToDo[head];
    verbose && printf(verbose, "Chunk #%d of %d (%d units)\n",
                                        count, nbrOfChunks, length(units));
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    unitsToDo <- unitsToDo[-head];

    # Read CDF structure
    verbose && enter(verbose, "Reading CDF list structure");
    srcUnits <- .readCdf(src, units=units, readGroupDirection=TRUE);
    verbose && exit(verbose);

    # Sanity check
    if (is.null(srcUnits)) {
      throw(sprintf("Failed to read %d units from CDF file.  This could be because you are running out of memory.  Try decreasing argument 'ram': %s", length(units), src));
    }

    # Sanity check
    if (length(srcUnits) != length(units)) {
      throw("Number of read CDF units does not equal number of requested units: ", length(srcUnits), " != ", length(units));
    }

    if (verbose && isVisible(verbose)) {
      cat(sprintf(" => RAM: %.fMB\n", object.size(srcUnits)/1024^2));
    }

    # Sanity check
    if (length(srcUnits) == 0) {
      throw("Internal error: While creating unique CDF, an empty list of CDF units was read.");
    }

    # Precalculate
    srcUnits <- lapply(srcUnits, FUN=function(unit) {
      groups <- .subset2(unit, "groups");
      groups <- lapply(groups, FUN=function(group) {
		nThisGroup <- length(.subset2(group,"pbase"));
        #group[fields] <- lapply(.subset(group, fields), FUN=.subset, fIdxs);
        #group$natoms <- nThisGroup
        #group$ncellsperatom <- as.integer(1);
        idxs <- idxOffset + seq_len(nThisGroup);
        idxs1 <- as.integer(idxs-1);  # Makes sure x & y are integers.
        y <-  idxs1 %/% ncols;
        group$y <- y;
        group$x <- idxs1 - ncols*y;
        idxOffset <<- idxOffset + nThisGroup;
        group;
      })
      #ncells <- length(groups)*nbrOfCellsPerField;
      unit$groups <- groups;
      #unit$ncells <- ncells;
      #unit$natoms <- ncells;
      #unit$ncellsperatom <- as.integer(1);
      unit;
    })

    #return(srcUnits)

    # Sanity check
    if (length(srcUnits) == 0) {
      throw("Internal error: While creating unique CDF, an empty list of CDF units is requested to be written.");
    }

#    verbose && str(verbose, srcUnits[1]);

    # Write regular units
    .writeCdfUnits(con=con, srcUnits, verbose=verbose2);
    # Not needed anymore
    srcUnits <- units <- NULL; # Not needed anymore

    count <- count + 1;
  } # while (length(unitsToDo) > 0)
  # Not needed anymore
  unitsToDo <- head <- fields <- count <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Close the connection
  close(con);
  con <- NULL;

  verbose && cat(verbose, "Temporary CDF file created:");
  verbose && cat(verbose, "Output pathname: ", pathnameT);
  verbose && print(verbose, file.info(pathnameT));

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Verifying the CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Verifying the written CDF");
  # Checking header
  header <- .readCdfHeader(pathnameT);

  # Validation array dimension
  if ((header$nrows != nrows) || (header$ncols != ncols)) {
    throw(sprintf("Failed to create a valid unique-cell CDF: The dimension of the written CDF does not match the intended one: (%d,%d) != (%d,%d)", header$nrows, header$ncols, nrows, ncols));
  }

  # Checking cell indices
  # Check CDF file in chunks
  nbrOfUnits <- header$probesets;
  chunkSize <- 10000*ram;
  nbrOfChunks <- ceiling(nbrOfUnits / chunkSize);
  for (kk in 1:nbrOfChunks) {
    verbose && printf(verbose, "Chunk %d of %d\n", kk, nbrOfChunks);
    from <- (kk-1)*chunkSize+1;
    to <- min(from+chunkSize, nbrOfUnits);
    cells <- .readCdfCellIndices(pathnameT, units=from:to);
    cells <- unlist(cells, use.names=FALSE);
    cells <- diff(cells);
    cells <- unique(cells);
    udcells <- as.integer(cells);
    if (!identical(udcells, 1:1)) {
      throw("Failed to create a valid unique-cell CDF: The cell indices are not contiguous: ", paste(udcells, collapse=", "));
    }
    # Not needed anymore
    cells <- udcells <- NULL;
  }
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=verbose);

  verbose && cat(verbose, "File pathname: ", pathname);
  verbose && print(verbose, file.info(pathname));
  verbose && exit(verbose);

  verbose && exit(verbose);

  # Return an AffymetrixCdfFile object for the new CDF
  cdfU <- newInstance(this, pathname);


  verbose && enter(verbose, "Final set of sanity checks");

  verbose && cat(verbose, "Number of units");
  stopifnot(nbrOfUnits(cdfU) == nbrOfUnits(this));

  verbose && cat(verbose, "Number of groups per unit");
  stopifnot(identical(nbrOfGroupsPerUnit(cdfU), nbrOfGroupsPerUnit(this)));

  verbose && cat(verbose, "Groups names per unit");
  stopifnot(identical(.readCdfGroupNames(getPathname(cdfU)), .readCdfGroupNames(getPathname(this))));

  verbose && cat(verbose, "Number of cells per unit group");
  stopifnot(identical(nbrOfCellsPerUnitGroup(cdfU), nbrOfCellsPerUnitGroup(this)));

  verbose && cat(verbose, "Consecutive ordering of cell indices");
  cells <- getCellIndices(cdfU, unlist=TRUE, useNames=FALSE);
  stopifnot(length(cells) <= nbrOfCells(cdfU));
  stopifnot(identical(unique(diff(cells)), 1L));

  ## Create checksum file
  cdfUZ <- getChecksumFile(cdfU)

  verbose && exit(verbose);

  cdfU;
}, private=TRUE) # createUniqueCdf()




setMethodS3("getUniqueCdf", "AffymetrixCdfFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Already a unique CDF?
  if (isUniqueCdf(this)) {
    return(this);
  }

  verbose && enter(verbose, "Retrieving unique CDF");

  # The full chiptype of the unique CDF
  chipType <- sprintf("%s,unique", getChipType(this));
  verbose && cat(verbose, "Unique chip type: ", chipType);

  # First, try to locate an existing unique CDF
  verbose && enter(verbose, "Locating unique CDF");
  pathname <- findByChipType(this, chipType=chipType, ...);
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && exit(verbose);

  if (is.null(pathname)) {
    verbose && enter(verbose, "Could not locate unique CDF. Will create one for chip type");
    res <- createUniqueCdf(this, ..., verbose=less(verbose));
    verbose && exit(verbose);
  } else {
    res <- byChipType(this, chipType=chipType, ...);
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE)


setMethodS3("isUniqueCdf", "AffymetrixCdfFile", function(this, ...) {
  hasTag(this, "unique");
})


# equals(getMainCdf(getUniqueCdf(cdf), cdf)) == TRUE
#setMethodS3("getMainCdf", "AffymetrixCdfFile", function(this, ...) {
#  # Argument 'this':
#  if (!isUniqueCdf(this)) {
#    throw("Cannot get the main CDF, because this CDF is not a unique CDF: ", getPathname(this));
#  }
#
#  chipType <- getChipType(this, fullname=FALSE);
#  tags <- getTags(this);
#  tags <- setdiff(tags, c("unique"));
#
#  # Try to locate the main CDF
#  cdf <- byChipType(this, chipType=chipType, tags=tags, ...);
#
#  cdf;
#}, protected=TRUE)



setMethodS3("getUnitGroupCellMapWithUnique", "AffymetrixCdfFile", function(this, units=NULL, ..., ugcMapM=NULL, verbose=FALSE) {
  # Argument 'this':
  if (isUniqueCdf(this)) {
    throw("Argument 'this' is already a unique CDF: ", getPathname(this));
  }

  # Argument 'units':
  units0 <- units;

  # Argument 'ugcMapM':
  if (is.null(ugcMapM)) {
    cdfM <- getUniqueCdf(this);
  } else if (!isUnitGroupCellMap(ugcMapM)) {
    throw("Argument 'ugcMapM' is not a UnitGroupCellMap");
  }

  # Merge maps with ordered and unique units (and the expand in the end)
  if (!is.null(units))
    units <- sort(unique(units));

  ugcMap <- list();

  # Get the (unit, group, cell) map for this CDF
  ugcMap[[1]] <- getUnitGroupCellMap(this, units=units, verbose=verbose);

  # Get the (unit, group, cell) map for the unique CDF
  if (is.null(ugcMapM)) {
    ugcMap[[2]] <- getUnitGroupCellMap(cdfM, units=units, verbose=verbose);
  } else {
    ugcMap[[2]] <- ugcMapM;
  }

  # Expand the unique CDF (unit, group, cell) map to the main one
  # (1) Create a (unit, group) -> hashcode map
  MAXGROUP <- 10000;
  ugHashcode <- lapply(ugcMap, FUN=function(map) {
    units  <- map[,"unit"];
    groups <- map[,"group"];
    hashcode <- MAXGROUP*units + groups;
    hashcode;
  })
  # (2) Map (unit, group)_unique to (unit, group)_main
  rr <- match(ugHashcode[[1]], ugHashcode[[2]]);
  # Not needed anymore
  ugHashcode <- NULL;
  mergedMap <- cbind(ugcMap[[1]], cellM=ugcMap[[2]][rr,"cell"]);
  # Not needed anymore
  ugcMap <- NULL;

  mergedMap;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2011-04-15 [HB]
# o Added more help to createUniqueCdf() for AffymetrixCdfFile.
# 2008-11-28 [HB]
# o BUG FIX: createUniqueCdf() used 'cdf' instead of 'this'.
# 2008-11-14 [MR]
# o created from Affymetrix.MONOCELL.R
############################################################################
