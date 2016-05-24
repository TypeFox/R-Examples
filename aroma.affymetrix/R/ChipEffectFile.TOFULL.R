setMethodS3("getExpandedCellMap", "ChipEffectFile", function(this, resetFields=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  cells <- getCellIndices(this, ..., verbose=verbose);
  names(cells) <- NULL;
  cells <- lapply(cells, FUN=unlist, use.names=FALSE);
  unitSizes <- sapply(cells, FUN=length);
  uUnitSizes <- unique(unitSizes);

  ceX <- clone(this);
  ## verbose && print(verbose, resetFields);
  for (field in resetFields) {
    ceX[[field]] <- FALSE;
  }

  cellsX <- getCellIndices(ceX, ..., verbose=verbose);
  names(cellsX) <- NULL;
  cellsX <- lapply(cellsX, FUN=unlist, use.names=FALSE);
  unitSizesX <- sapply(cellsX, FUN=length);

  resizeFactors <- unitSizesX / unitSizes;
  verbose && print(verbose, table(resizeFactors));

  map <- rep(as.integer(NA), times=nbrOfCells(this));

  verbose && cat(verbose, "Unit sizes:");
  verbose && print(verbose, table(unitSizes));
#  # Not needed anymore
#  unitSizesX <- unitSizes <- NULL;

  for (size in uUnitSizes) {
    verbose && enter(verbose, "Unit size: ", size);
    units <- which(unitSizes == size);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);
    times <- resizeFactors[units];
    verbose && cat(verbose, "Times:");
    verbose && str(verbose, times);
    units1 <- rep(units, times=times);
    cells1 <- unlist(cells[units1], use.names=FALSE);
    verbose && cat(verbose, "Cells 1:");
    verbose && str(verbose, cells1);
    cells2 <- unlist(cellsX[units], use.names=FALSE);
    verbose && cat(verbose, "Cells 2:");
    verbose && str(verbose, cells2);
    map[cells2] <- cells1;
    # Not needed anymore
    cells1 <- cells2 <- units <- NULL;
    verbose && exit(verbose);
  }
  # Not needed anymore
  cells <- cellsX <- NULL;

  stopifnot(length(map) == nbrOfCells(this));

  map;
}, protected=TRUE)


setMethodS3("getExpandedCellMap", "SnpChipEffectFile", function(this, resetFields=NULL, ...) {
  NextMethod("getExpandedCellMap", resetFields=c("mergeStrands", resetFields));
}, protected=TRUE)


setMethodS3("getExpandedCellMap", "CnChipEffectFile", function(this, resetFields=NULL, ...) {
  NextMethod("getExpandedCellMap", resetFields=c("combineAlleles", resetFields));
}, protected=TRUE)


setMethodS3("getCellMapForMainCdf", "ChipEffectFile", function(this, ..., verbose=FALSE) {
  cdfM <- getCdf(this);
  cdf <- getMainCdf(cdfM);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a map from the indices of the monocell to the main CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  map <- getUnitGroupCellMapWithMonocell(cdf, verbose=verbose);
  readMap <- map[,"cellM"];
  o <- order(readMap);
  writeMap <- map[o,"cell"];
  # Not needed anymore
  map <- NULL;
  readMap <- readMap[o];
  # Not needed anymore
  o <- NULL;

  map2 <- getExpandedCellMap(this, verbose=verbose);
  readMap <- map2[readMap];

  readMap2 <- rep(as.integer(NA), length=nbrOfCells(cdf));
  readMap2[writeMap] <- readMap;
  # Not needed anymore
  readMap <- writeMap <- NULL;

  # Sanity check
  stopifnot(length(readMap2) == nbrOfCells(cdf));

  readMap2;
}, protected=TRUE)


setMethodS3("writeAsFullCelFile", "ChipEffectFile", function(this, name=getName(this), tags="*", ..., cells=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (!is.null(tags)) {
    tags <- Arguments$getTags(tags, collapse=NULL);
    idx <- which(tags == "*");
    if (length(idx) > 0)
      tags <- insert(tags[-idx], at=idx, getTags(this));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create empty CEL file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdfM <- getCdf(this);
  cdf <- getMainCdf(cdfM);

  # Allocate a temporary CEL file (note suffix argument)
  cf <- AffymetrixCelFile$allocateFromCdf(cdf, name=name, tags=tags, suffix=".CEL.tmp", ..., verbose=less(verbose, 10));
  verbose && print(verbose, cf);
  pathnameT <- getPathname(cf);

  # Get the index map that maps monocell-CDF cells to main-CDF cells
  if (is.null(cells)) {
    cells <- getCellMapForMainCdf(this, verbose=less(verbose, 10));
    verbose && str(verbose, cells);
  }

  # Read data
  fields <- c("intensities", "stdvs", "pixels");
  data <- readRawData(this, indices=cells, fields=fields, verbose=less(verbose, 10));
  verbose && str(verbose, data);

  # Identify cells that have data
  mainCells <- which(is.finite(cells));
  data <- data[mainCells,,drop=FALSE];
  verbose && str(verbose, data);

  .updateCel(pathnameT, indices=mainCells, intensities=data);

  # Rename from temporary to final filename (see above)
  pathname <- popTemporaryFile(pathnameT, verbose=verbose);

  res <- AffymetrixCelFile(pathname);

  attr(res, "cells") <- cells;

  res;
}, protected=TRUE);



setMethodS3("getAsFullCelFile", "ChipEffectFile", function(this, name=getName(this), tags="*", path=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (is.null(path)) {
    rootPath <- "probeData";
    path <- getPath(this);
    chipType <- basename(path);
    dataSet <- basename(getParent(path));
    path <- file.path(rootPath, dataSet, chipType);
    path <- Arguments$getWritablePath(path);
  }
  path <- Arguments$getWritablePath(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (!is.null(tags)) {
    tags <- Arguments$getTags(tags, collapse=NULL);
    idx <- which(tags == "*");
    if (length(idx) > 0) {
      tags <- insert(tags[-idx], at=idx, getTags(this));
    }
    tags <- tags[(tags != "chipEffects")];
  }


  fullname <- paste(c(name, tags), collapse=",");
  filename <- sprintf("%s.CEL", fullname);
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  if (isFile(pathname)) {
    res <- AffymetrixCelFile(pathname);
  } else {
    res <- writeAsFullCelFile(this, name=name, tags=tags, path=path, ...,
                                                 verbose=less(verbose, 2));
  }

  res;
}, protected=TRUE);


############################################################################
# HISTORY:
# 2008-03-31
# o Now getAsFullCelFile() writes to probeData/.
# 2008-03-18
# o Added getAsFullCelFile() and writeAsFullCelFile() for ChipEffectFile.
############################################################################
