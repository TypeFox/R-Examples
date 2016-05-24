###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod createMonocellCdf
# @aliasmethod createMonoCell
#
# @title "Creates a mono-cell version of the CDF"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{chipType}{The chip type of the new CDF.}
#   \item{tags}{Tags added to the chip type of the new CDF.}
#   \item{path}{The path where to store the new CDF file.}
#   \item{nbrOfCellsPerField}{Number of cells per group field the new CDF
#      should have.}
#   \item{...}{Additional arguments passed to @see "affxparser::writeCdf".}
#   \item{ram}{A @double saying if more or less units should be converted
#      per chunk.  A smaller value uses less memory.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the monocell CDF as an @see "AffymetrixCdfFile" object.
# }
#
# @author "HB, KS"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("createMonocellCdf", "AffymetrixCdfFile", function(this, chipType=getChipType(this), tags="monocell", path=NULL, nbrOfCellsPerField=1, ..., ram=NULL, overwrite=FALSE, verbose=TRUE) {
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
          # nindices <- length(group$indices);
          nindices <- length(.subset2(group, "indices"));
          head <- 1:nindices;
          # idxs <- cells[head];
          idxs <- .subset(cells, head);
          # cells <- cells[<tail>];
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
        # group <- units[[kk]]
        group <- .subset2(units, kk);
        # Number of cells in this group
        # nindices <- length(group$indices);
        nindices <- length(.subset2(group, "indices"));
        head <- 1:nindices;
        # idxs <- cells[head];
        idxs <- .subset(cells, head);
        # cells <- cells[<tail>];
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
  # Argument 'nbrOfCellsPerField':
  nbrOfCellsPerField <- Arguments$getIndices(nbrOfCellsPerField);

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  verbose2 <- as.integer(verbose)-5;  # For 'affxparser' calls.


  # Already a monocell CDF?
  if (isMonocellCdf(this)) {
    return(this);
  }

  verbose && enter(verbose, "Creating monocell CDF");

  verbose && cat(verbose, "Chip type: ", getChipType(this));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate source CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Validate (main) CDF");
  validate(this);
  verbose && exit(verbose);


  # Get the pathname of the source
  src <- getPathname(this);
  src <- Arguments$getReadablePathname(src);

  # Create the pathname of the destination CDF
  if (is.null(path)) {
    mainChipType <- gsub("[,].*", "", chipType);
    path <- filePath("annotationData", "chipTypes", mainChipType);
  }

  name <- paste(c(chipType, tags), collapse=",");
  filename <- sprintf("%s.CDF", name);
  pathname <- Arguments$getWritablePathname(filename, path=path);
  pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=!overwrite);

  if (overwrite && isFile(pathname)) file.remove(pathname)

  # Write to a temporary file first, and rename it when we know it's complete
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

  # Assure source and destination is not identical
  if (identical(src, pathnameT)) {
    throw("Cannot create CDF file. Destination is same as source: ", src);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fields to be kept
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Number of cells to keep in each group field
  fIdxs <- 1:nbrOfCellsPerField;

  verbose && printf(verbose, "Number of cells per group field: %d\n",
                                                       nbrOfCellsPerField);


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

  # Get the number of cells per unit
  nbrOfCellsPerUnit <- nbrOfGroupsPerUnit * nbrOfCellsPerField;
  verbose && cat(verbose, "Number of cells per unit:");
  verbose && summary(verbose, nbrOfCellsPerUnit);

  # Total number of cells
  nbrOfCells <- sum(nbrOfCellsPerUnit);

  # Number of units
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

    unitsToDo <- unitsToDo[-head];

    # Read CDF structure
#    srcUnits <- readCdf(src, units=units);

# readGroupDirection is needed in order for writeCdf() to work. /KS 18/12/06
    verbose && enter(verbose, "Reading CDF list structure");
    srcUnits <- .readCdf(src, units=units, readGroupDirection=TRUE);
    verbose && exit(verbose);
    if (is.null(srcUnits)) {
      throw(sprintf("Failed to read %d units from CDF file.  This could be because you are running out of memory.  Try decreasing argument 'ram': %s", length(units), src));
    }

    if (length(srcUnits) != length(units)) {
      throw("Number of read CDF units does not equal number of requested units: ", length(srcUnits), " != ", length(units));
    }

    if (verbose && isVisible(verbose)) {
      cat(sprintf(" => RAM: %.fMB\n", object.size(srcUnits)/1024^2));
    }

    if (length(srcUnits) == 0) {
      throw("Internal error: While creating monocell, an empty list of CDF units was read.");
    }

    # 15/12/06: readCdf() no longer returns a list
    # with element "blocks"; instead, it returns a list
    # with element "groups". /KS

    # Precalculate
    srcUnits <- lapply(srcUnits, FUN=function(unit) {
      groups <- .subset2(unit, "groups");
      groups <- lapply(groups, FUN=function(group) {
        group[fields] <- lapply(.subset(group, fields), FUN=.subset, fIdxs);
        group$natoms <- nbrOfCellsPerField;
        group$ncellsperatom <- as.integer(1);
        idxs <- idxOffset + fIdxs;
        idxs1 <- as.integer(idxs-1);  # Makes sure x & y are integers.
        y <-  idxs1 %/% ncols;
        group$y <- y;
        group$x <- idxs1 - ncols*y;
#        group$indices <- idxs; # Not needed (saves ~225Mb for 250K Nsp!)
        idxOffset <<- idxOffset + nbrOfCellsPerField;
        group;
      })
      ncells <- length(groups)*nbrOfCellsPerField;
      unit$groups <- groups;
      unit$ncells <- ncells;
      unit$natoms <- ncells;
      unit$ncellsperatom <- as.integer(1);
      unit;
    })

    if (length(srcUnits) == 0) {
      throw("Internal error: While creating monocell, an empty list of CDF units is requested to be written.");
    }

#    verbose && str(verbose, srcUnits[1]);

    # Write regular units
    .writeCdfUnits(con=con, srcUnits, verbose=verbose2);
    # Not needed anymore
    srcUnits <- units <- NULL; # Not needed anymore

    count <- count + 1;
  } # while (length(unitsToDo) > 0)
  # Not needed anymore
  unitsToDo <- head <- fields <- fIdxs <- count <- NULL;

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
  if ((header$nrows != nrows) || (header$ncols != ncols)) {
    throw(sprintf("Failed to create a valid mono-cell CDF: The dimension of the written CDF does not match the intended one: (%d,%d) != (%d,%d)", header$nrows, header$ncols, nrows, ncols));
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
      throw("Failed to create a valid mono-cell CDF: The cell indices are not contiguous: ", paste(udcells, collapse=", "));
    }
    # Not needed anymore
    cells <- udcells <- NULL;
  }
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Renaming temporary file
  popTemporaryFile(pathnameT, verbose=verbose);

  verbose && cat(verbose, "File pathname: ", pathname);
  verbose && print(verbose, file.info(pathname));
  verbose && exit(verbose);

  verbose && exit(verbose);

  # Return an AffymetrixCdfFile object for the new CDF
  cdfM <- newInstance(this, pathname)

  ## Create checksum file
  cdfMZ <- getChecksumFile(cdfM)

  cdfM
}, private=TRUE) # createMonocellCdf()


setMethodS3("getMonocellCdf", "AffymetrixCdfFile", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Already a monocell CDF?
  if (isMonocellCdf(this)) {
    return(this);
  }

  verbose && enter(verbose, "Retrieving monocell CDF");

  # The full chiptype of the monocell CDF
  chipType <- sprintf("%s,monocell", getChipType(this));
  verbose && cat(verbose, "Monocell chip type: ", chipType);

  # First, try to locate an existing monocell CDF
  verbose && enter(verbose, "Locating monocell CDF");
  pathname <- findByChipType(this, chipType=chipType, ..., verbose=less(verbose,20));
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && exit(verbose);

  if (force || is.null(pathname)) {
    verbose && enter(verbose, "Could not locate monocell CDF. Will create one for chip type");
    res <- createMonocellCdf(this, ..., overwrite=force, verbose=less(verbose));
    verbose && exit(verbose);
  } else {
    res <- byChipType(this, chipType=chipType, ...);
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE)


setMethodS3("isMonocellCdf", "AffymetrixCdfFile", function(this, ...) {
  hasTag(this, "monocell");
})



# equals(getMainCdf(getMonocellCdf(cdf), cdf)) == TRUE
setMethodS3("getMainCdf", "AffymetrixCdfFile", function(this, ...) {
  # Argument 'this':
  if (!isMonocellCdf(this) & !isUniqueCdf(this)) {
    throw("Cannot get the main CDF, because this CDF is not a monocell/unique CDF: ", getPathname(this));
  }

  chipType <- getChipType(this, fullname=FALSE);
  tags <- getTags(this);
  tags <- setdiff(tags, c("monocell"));
  tags <- setdiff(tags, c("unique"));

  # Try to locate the main CDF
  cdf <- byChipType(this, chipType=chipType, tags=tags, ...);

  cdf;
}, protected=TRUE)



setMethodS3("getUnitGroupCellMapWithMonocell", "AffymetrixCdfFile", function(this, units=NULL, ..., ugcMapM=NULL, verbose=FALSE) {
  # Argument 'this':
  if (isMonocellCdf(this)) {
    throw("Argument 'this' is already a monocell CDF: ", getPathname(this));
  }

  # Argument 'units':
  units0 <- units;

  # Argument 'ugcMapM':
  if (is.null(ugcMapM)) {
    cdfM <- getMonocellCdf(this);
  } else if (!isUnitGroupCellMap(ugcMapM)) {
    throw("Argument 'ugcMapM' is not a UnitGroupCellMap");
  }

  # Merge maps with ordered and unique units (and the expand in the end)
  if (!is.null(units))
    units <- sort(unique(units));

  ugcMap <- list();

  # Get the (unit, group, cell) map for this CDF
  ugcMap[[1]] <- getUnitGroupCellMap(this, units=units, verbose=verbose);

  # Get the (unit, group, cell) map for the monocell CDF
  if (is.null(ugcMapM)) {
    ugcMap[[2]] <- getUnitGroupCellMap(cdfM, units=units, verbose=verbose);
  } else {
    ugcMap[[2]] <- ugcMapM;
  }

  # Expand the monocell CDF (unit, group, cell) map to the main one
  # (1) Create a (unit, group) -> hashcode map
  MAXGROUP <- 10000;
  ugHashcode <- lapply(ugcMap, FUN=function(map) {
    units  <- map[,"unit"];
    groups <- map[,"group"];
    hashcode <- MAXGROUP*units + groups;
    hashcode;
  })
  # (2) Map (unit, group)_monocell to (unit, group)_main
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
# 2015-01-23
# o Added argument 'force' to getMonocellCdf() for AffymetrixCdfFile.
#   If @TRUE and an monocell CDF already exists, then it will be re-created
#   and overwritten, otherwise an error will be thrown.
# 2012-10-18
# o ROBUSTNESS: Now createMonocellCdf() for AffymetrixCdfFile validates
#   the CDF before trying to create the monocell CDF.  This should catch
#   for instance invalid Brainarray custom CDFs, cf. ...
# 2011-04-15
# o CLEANUP: Dropped argument 'sep' of createMonocellCdf() and
#   createUniqueCdf() for AffymetrixCdfFile.
# 2011-01-09
# o BUG FIX: createMonocellCdf(..., verbose=FALSE) for AffymetrixCdfFile
#   would still create verbose output.
# 2008-12-15 [MR]
# o Now getMainCdf() works for "unique" (as well as "monocell") CDFs.
# 2008-07-22
# o MEMORY OPTIMIZATION: Now the validation part of createMonocellCdf() for
#   AffymetrixCdfFile is also sensitive to the 'ram' argument.  The
#   validation is now also down in smaller chunks.
# 2008-03-18
# o Added getUnitGroupCellMapWithMonocell() to AffymetrixCdfFile.
# o Added getMainCdf() for AffymetrixCdfFile to get main CDF from a
#   monocell CDF.
# 2008-02-20
# o Now getMonocellCdf() and createMonocellCdf() quietly return the input
#   CDF if it already is a monocell CDF.
# o Added isMonocellCdf().
# o Renamed getMonoCell() to getMonocellCdf() and createMonoCell() to
#   createMonocellCdf(), because the former had strange names.
# 2007-12-11
# o Now getMonoCell() will create the monocell CDF, if missing.
# 2007-07-13
# o Now the createMonoCell() first writes a temporary CDF named *.cdf.tmp,
#   and then rename it to *.cdf when it has been validated.
# 2007-03-08
# o Added getMonoCell().
# 2007-02-21 /HB + KS
# o BUG FIX: When creating a monocell, the output did not strip of the tags
#   from the chip type, e.g. annotationData/Foo,core/For,core,monocell.cdf
#   instead of annotationData/Foo/For,core,monocell.cdf.
# 2007-02-08
# o Now monocell CDF are names <chipType>,monocell.cdf.  Before a dash was
#   used instead of a comma. This new style is more in line with the
#   <name>,<tags> naming convention used elsewhere in the package.
# 2007-02-06
# o Now monocell CDF are stored under annotationData/chipTypes/<chipType>/.
# 2007-01-10
# o Now createMonoCell() create the CDF in chunks, that is, in constant
#   memory.
# o Reordered internally in createMonoCell() preparing for code to read
#   *and* write monocell CDFs in chunks.  It should not be too hard.
#   We need to update affxparser with writeCdfHeader(), writeCdfQcUnits()
#   and writeCdfUnits(), and are basically already in there, but as
#   private functions.
# o Removed some unnecessary group fields in 'destUnits' saving us approx
#   220-230Mb out of 1.9GB for the Mapping250K_Nsp.  Still need to find
#   another solution to get down below 1GB. One thing that takes up a lot
#   of memory is that the unit and group directions are stored as character
#   strings and not integers.
# o BUG FIX: The most recent createMonoCell() would create CDFs with all
#   cell indices being equal to one.  Added more verbose output and some
#   garbage collection to this function too.
# 2007-01-06
# o Added argument 'force' to getCellIndices().
# o Now getCellIndices() only caches object < 10 MB RAM.
# o Optimized identifyCells() to only cache data in rare cases.
# 2006-12-18 /KS
# o Made global replacement "block" -> "group".
# 2006-12-14
# o Added convertUnits().
# 2006-09-27
# o Now fromFile() tries to create an instance of the subclasses (bottom up)
#   first.  This will make it possible to automatically define SNP CDFs.
# 2006-09-26
# o Now getGenomeInformation() and getSnpInformation() ignores suffices of
#   the chip-type string. This makes it possible to retrive annotation data
#   also for custom chips.
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-16
# o Added getGenomeInformation() and stextChipType().
# 2006-09-14
# o BUG FIX: Fractional value of 'indices' of identifyCells().
# 2006-09-10
# o BUG FIX: createMonoCell() where resetting the cell indices for each
#   chunk.
# o Simple benchmarking of createMonoCell(): IBM Thinkpad A31 1.8GHz 1GB:
#   Mapping50K_Hind to mono cells CDF takes ~13 mins.  Again, it is
#   writeCdf() that is slow.  KH is working on improving this.
# 2006-09-08
# o Added equals() to compare to CDF object.
# o Added convert() to convert a CDF into another version by convertCdf().
# 2006-08-25
# o Added getFirstCellIndices().  This may be used by methods to identify
#   unit or unit groups for which no probe signals have been assigned yet.
# 2006-08-24
# o Added reconstruct().
# o Added Rdoc comments.
# 2006-08-19
# o Added getUnitSizes().  Useful in order to store parameter estimates
#   of probeset-summary models.
# 2006-08-11
# o Created.
############################################################################
