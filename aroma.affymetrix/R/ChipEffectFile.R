###########################################################################/**
# @RdocClass ChipEffectFile
#
# @title "The ChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ParameterCelFile".}
#   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getChipEffectSet()} method for the @see "ProbeLevelModel" class.
#   An object of this class is typically part of a @see "ChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("ChipEffectFile", function(..., probeModel=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  this <- extend(ParameterCelFile(...), "ChipEffectFile",
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )

  setEncodeFunction(this, function(groupData, ...) {
    theta <- .subset2(groupData, "theta");
    stdvs <- .subset2(groupData, "sdTheta");
    outliers <- .subset2(groupData, "thetaOutliers");
    pixels <- NULL;
    if (!is.null(outliers))
      pixels <- -as.integer(outliers);

    res <- list();
    if (!is.null(theta))
      res$intensities <- theta;
    if (!is.null(stdvs))
      res$stdvs <- stdvs;
    if (!is.null(pixels))
      res$pixels <- pixels;

    res;
  })

  ## Not safe, at least not what I know of. /HB 2007-08-17
##   setEncodeFunction(this, function(groupData, ...) {
##     groupData[[3]] <- -as.integer(.subset2(groupData, 3));
##     names(groupData) <- c("intensities", "stdvs", "pixels");
##     groupData;
##   })


  setDecodeFunction(this, function(groupData, ...) {
    res <- list();
    if (!is.null(groupData$intensities))
      res$theta <- groupData$intensities;
    if (!is.null(groupData$stdvs))
      res$sdTheta <- groupData$stdvs;
    if (!is.null(groupData$pixels))
      res$thetaOutliers <- as.logical(-groupData$pixels);
    res;
  })

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("as.character", "ChipEffectFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  s;
}, protected=TRUE)


setMethodS3("getParameters", "ChipEffectFile", function(this, ...) {
  params <- NextMethod("getParameters");
  params$probeModel <- this$probeModel;
  params;
}, protected=TRUE)


setMethodS3("createParamCdf", "ChipEffectFile", function(static, sourceCdf, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Creating CDF for chip effects");
  verbose && cat(verbose, "Source chip type: ", getChipType(sourceCdf));
  verbose && cat(verbose, "Source CDF: ", getPathname(sourceCdf));

  # Search for existing monocell CDF
  for (sep in c(",", "-")) {
    chipType <- paste(getChipType(sourceCdf), "monocell", sep=sep);
    verbose && cat(verbose, "Looking for chip type: ", chipType);
    pathname <- AffymetrixCdfFile$findByChipType(chipType);
    if (!is.null(pathname)) {
      verbose && cat(verbose, "Found: ", pathname);
      break;
    }
  }
  # Warn about deprecated filname <chipType>-monocell.
  if (!is.null(pathname) && (sep == "-")) {
    msg <- paste("Deprecated filename of monocell CDF detected (uses dash instead of comma): ", pathname);
    verbose && cat(verbose, msg);
    throw(msg);
  }

  if (is.null(pathname)) {
    verbose && cat(verbose, "Pathname: Not found!");
    verbose && cat(verbose, "Will create CDF for the chip-effect files from the original CDF. NOTE: This will take several minutes or more!");
    verbose && enter(verbose, "Creating CDF");
    cdf <- createMonocellCdf(sourceCdf, verbose=less(verbose));
    verbose && exit(verbose);
  } else {
    verbose && cat(verbose, "Pathname: ", pathname);
    cdf <- AffymetrixCdfFile$fromFile(pathname);
  }
  verbose && exit(verbose);

  cdf;
}, static=TRUE, private=TRUE)




setMethodS3("readUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, ..., force=FALSE, cache=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # Check for cached data
  key <- list(method="readUnits", class=class(this)[1],
              pathname=getPathname(this),
              cdf=cdf, units=units, ...);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="readUnits", pathname=getPathname(this), cdf=cdf, units=units, ...);
  }
  id <- getChecksum(key);
  res <- this$.readUnitsCache[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.ChipEffectFile(): Returning cached data");
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve the data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(cdf)) {
    cdf <- getCellIndices(this, units=units, verbose=less(verbose));
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  res <- NextMethod("readUnits", cdf=cdf, force=force, verbose=less(verbose));

  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.ChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
  }

  res;
})



###########################################################################/**
# @RdocMethod getCellIndices
#
# @title "Retrieves tree list of cell indices for a set of units"
#
# \description{
#   @get "title" from the associated CDF.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to \code{getCellIndices()}
#             of @see "AffymetrixCdfFile".}
#  \item{.cache}{Ignored.}
# }
#
# \value{
#   Returns a @list structure, where each element corresponds to a unit.
#   If argument \code{unlist=TRUE} is passed, an @integer @vector is returned.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getCellIndices", "ChipEffectFile", function(this, ..., .cache=TRUE) {
  cdf <- getCdf(this);
  getCellIndices(cdf, ...);
}, protected=TRUE)



setMethodS3("updateUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, data, ...) {
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  NextMethod("updateUnits", cdf=cdf, data=data);
}, private=TRUE);



setMethodS3("findUnitsTodo", "ChipEffectFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Identifying non-fitted units in chip-effect file");

  verbose && cat(verbose, "Pathname: ", getPathname(this));


  idxs <- NULL;
  if (is.null(units)) {
    # Look up chip-type and parameter specific but data set independent data
    cdf <- getCdf(this);
    chipType <- getChipType(cdf);
    key <- list(method="findUnitsTodo", class=class(this)[1],
                chipType=chipType, params=getParameters(this));
    dirs <- c("aroma.affymetrix", chipType);
    if (!force) {
      idxs <- loadCache(key, dirs=dirs);
      if (!is.null(idxs))
        verbose && cat(verbose, "Found indices cached on file");
    }
  }

  if (is.null(idxs)) {
    verbose && enter(verbose, "Identifying CDF units");

    units0 <- units;
    if (is.null(units)) {
      cdf <- getCdf(this);
      units <- seq_len(nbrOfUnits(cdf));
    }
    nbrOfUnits <- length(units);

    idxs <- lapplyInChunks(units, function(unitsChunk, ...) {
      verbose && enter(verbose, "Reading CDF cell indices");
      idxsChunk <- getCellIndices(this, units=unitsChunk, force=TRUE, verbose=less(verbose));
      # Save memory 100% -> 94%
      names(idxsChunk) <- NULL;
      verbose && exit(verbose);

      verbose && enter(verbose, "Extracting first CDF group for each unit");
      # Save memory 100% -> 94% -> 5% (all the nested named lists costs!)
      idxsChunk <- lapply(idxsChunk, FUN=function(unit) {
        groups <- .subset2(unit, "groups");
        fields <- .subset2(groups, 1);
        .subset2(fields, 1);
      })
      verbose && exit(verbose);

      gc <- gc();

      idxsChunk;
    }, chunkSize=100e3, useNames=FALSE, verbose=verbose);

    # Not needed anymore
    units <- units0;
    # Not needed anymore
    units0 <- NULL;

    idxs <- unlist(idxs, use.names=FALSE);
    gc <- gc();
    verbose && print(verbose, gc);

    # Verify correctness
    if (length(idxs) != nbrOfUnits) {
      throw("Internal error: Expected ", nbrOfUnits, " cell indices, but got ", length(idxs), ".");
    }

    if (is.null(units)) {
      verbose && enter(verbose, "Saving to file cache");
      saveCache(idxs, key=key, dirs=dirs);
      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  }


  # Read one cell from each unit
  verbose && enter(verbose, "Reading data for these ", length(idxs), " cells");
  value <- .readCel(getPathname(this), indices=idxs, readIntensities=FALSE,
                   readStdvs=TRUE, readPixels=FALSE)$stdvs;
  verbose && exit(verbose);


  # Identify units for which the stdvs <= 0.
  value <- which(value <= 0);
  if (!is.null(units))
    value <- units[value];
  verbose && cat(verbose, "Looking for stdvs <= 0 indicating non-estimated units:");
  verbose && str(verbose, value);

  verbose && exit(verbose);

  value;
})



###########################################################################/**
# @RdocMethod getUnitGroupCellMap
# @aliasmethod getCellMap
#
# @title "Gets a (unit, group, cell) index map"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units for which the map should be returned.
#      If @NULL, all units are considered.}
#   \item{force}{If @TRUE, cached cell indices are ignored.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @data.frame with @integer columns \code{unit}, \code{group},
#  and \code{cell}.
# }
#
# \examples{\dontrun{
#      unit group cell
#    1  104     1  335
#    2  104     2  336
#    3  105     1  337
#    4  105     2  338
#    5  105     3  339
#    6  105     4  340
# }}
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getUnitGroupCellMap", "ChipEffectFile", function(this, units=NULL, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (inherits(units, "UnitGroupCellMap")) {
    return(units);
  } else if (is.null(units)) {
  } else if (is.list(units)) {
  } else {
    units <- Arguments$getIndices(units);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving unit-to-cell map");

  # Special case: requesting zero units?
  if (length(units) == 0 && !is.null(units)) {
    map <- data.frame(unit=integer(0), group=integer(0), cell=integer(0));
    class(map) <- c("UnitGroupCellMap", class(map));
    verbose && exit(verbose);
    return(map);
  }


  # Get the CDF
  cdf <- getCdf(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  useFileCache <- (is.null(units) || (!is.list(units) && length(units) > 10000));
  useFileCache <- (is.null(units) || (!is.list(units) && length(units) > 100));
  if (useFileCache) {
    chipType <- getChipType(cdf);
    # Look up chip-type and parameter specific but data set independent data
    key <- list(method="getUnitGroupCellMap", class=class(this)[1],
                chipType=chipType, params=getParameters(this),
                units=units);
    dirs <- c("aroma.affymetrix", chipType);
    if (!force) {
      map <- loadCache(key, dirs=dirs);
      if (!is.null(map)) {
        verbose && cat(verbose, "Found (unit,group,cell) map cached on file");
        verbose && exit(verbose);
        return(map);
      }
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Is 'units' already a CDF list?
  if (is.list(units)) {
    # No fancy validation for now.
    cells <- units;
    units <- indexOf(cdf, names=names(units));
    if (any(is.na(units))) {
      throw("Argument 'units' is of unknown structure.");
    }
    verbose && enter(verbose, "Argument 'cells' is already a CDF cell-index structure");

    # Get the unit names
    unitNames <- names(cells);

    # Get the number of groups per unit
    unitSizes <- lapply(cells, FUN=function(unit) {
      length(.subset2(unit, "groups"));
    });
    unitSizes <- unlist(unitSizes, use.names=FALSE);

    cells <- unlist(cells, use.names=FALSE);
  } else {
    verbose && enter(verbose, "Retrieving cell indices for specified units");
    if (is.null(units))
      units <- seq_len(nbrOfUnits(cdf));
    chunks <- splitInChunks(units, chunkSize=100e3);
    nbrOfChunks <- length(chunks);
    nbrOfUnits <- length(units);

    # Store unit names and unit sizes
    unitNames <- vector("character", nbrOfUnits);
    unitSizes <- vector("integer", nbrOfUnits);

    # Store cell indices first chunk-by-chunk, then as a vector.
    cells <- vector("list", nbrOfChunks);

    offset <- 0;
    for (kk in seq_len(nbrOfChunks)) {
      verbose && printf(verbose, "Chunk #%d of %d\n", kk, length(chunks));
      chunk <- chunks[[kk]];
      chunks[[kk]] <- NA;
      cells0 <- getCellIndices(this, units=chunk, force=force, .cache=FALSE, verbose=less(verbose));
      idxs <- offset + seq_len(length(chunk));
      offset <- offset + length(chunk);
      # Not needed anymore
      chunk <- NULL;
      unitNames[idxs] <- names(cells0);
      names(cells0) <- NULL;
      unitSizes0 <- lapply(cells0, FUN=function(unit) {
        length(.subset2(unit, "groups"));
      });
      unitSizes[idxs] <- unlist(unitSizes0, use.names=FALSE);
      # Not needed anymore
      unitSizes0 <- NULL;
      cells[[kk]] <- unlist(cells0, use.names=FALSE);
      # Not needed anymore
      cells0 <- idxs <- NULL;
    }
    # Not needed anymore
    chunks <- NULL;
    gc <- gc();
    verbose && print(verbose, gc);
  }

  cells <- unlist(cells, use.names=FALSE);
  gc <- gc();
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating return data frame");
  uUnitSizes <- sort(unique(unitSizes));
  verbose && cat(verbose, "Unique number of groups per unit: ",
                                        paste(uUnitSizes, collapse=","));
  verbose && cat(verbose, "Number of units: ", length(unitNames));

  if (is.null(units))
    units <- seq_len(nbrOfUnits(cdf));

  # The following is too slow:
  #  groups <- sapply(unitSizes, FUN=function(n) seq_len(n));

  # Instead, updated size by size
  verbose && printf(verbose, "Allocating matrix of size %dx%d.\n",
                                     max(uUnitSizes), length(unitNames));
  naValue <- as.integer(NA);
  units2 <- groups <- matrix(naValue, nrow=max(uUnitSizes), ncol=length(unitNames));
  for (size in uUnitSizes) {
    # Identify units with a certain number of groups
    cc <- which(unitSizes == size);
    # Assign group indices to the groups
    seq <- seq_len(size);
    groups[seq,cc] <- seq;
    # Assign unit indices to ditto
    units2[seq,cc] <- rep(units[cc], each=size);
    # Not needed anymore
    seq <- NULL;
    gc <- gc();
  }
  keep <- !is.na(groups);
  groups <- groups[keep];
  units2 <- units2[keep];
  # Not needed anymore
  keep <- NULL;
  gc <- gc();

  map <- data.frame(unit=units2, group=groups, cell=cells);
  verbose && exit(verbose);

  verbose && exit(verbose);

  class(map) <- c("UnitGroupCellMap", class(map));

  if (useFileCache) {
    verbose && enter(verbose, "Saving to file cache");
    saveCache(map, key=key, dirs=dirs);
    verbose && exit(verbose);
  }

  map;
}, private=TRUE)



setMethodS3("getUnitGroupCellChromosomePositionMap", "ChipEffectFile", function(this, units=NULL, chromosomes=NULL, orderByPosition=TRUE, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':
  ugcMap <- NULL;
  if (is.null(units)) {
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- ugcMap[,"unit"];
  }
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
  }
  units0 <- units;

  # Get the genome position information
  gi <- getGenomeInformation(cdf);

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
  chipType <- getChipType(cdf);
  key <- list(method="getUnitGroupCellChromosomePositionMap",
              class=class(this)[1],
              chipType=chipType, units=units, ugcMap=ugcMap,
              chromosomes=chromosomes, orderByPosition=orderByPosition);
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
  if (!isUnitGroupCellMap(ugcMap)) {
    ugcMap <- getUnitGroupCellMap(this, units=units, force=force, verbose=less(verbose, 10));
    verbose && cat(verbose, "(unit, group, cell) map:");
    verbose && str(verbose, ugcMap);
  }

  # Get the (chromosome, position) map
  cpMap <- getData(gi, units=ugcMap[,"unit"], force=force, verbose=less(verbose, 10));
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
}, private=TRUE)




setMethodS3("getDataFlat", "ChipEffectFile", function(this, units=NULL, fields=c("theta", "sdTheta", "outliers"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving data as a flat data frame");

  # Get unit-to-cell map
  suppressWarnings({
    map <- getUnitGroupCellMap(this, units=units, ..., verbose=less(verbose));
  })

  verbose && enter(verbose, "Reading data fields");
  celFields <- c(theta="intensities", sdTheta="stdvs", outliers="pixels");
  suppressWarnings({
    data <- getData(this, indices=map[,"cell"], fields=celFields[fields]);
  })
  rownames(data) <- seq_len(nrow(data));  # Work around?!? /HB 2006-11-28

  # Decode
  names <- colnames(data);
  names <- gsub("intensities", "theta", names);
  names <- gsub("stdvs", "sdTheta", names);
  names <- gsub("pixels", "outliers", names);
  colnames(data) <- names;
  verbose && str(verbose, data);
  if ("outliers" %in% names) {
    data[,"outliers"] <- as.logical(-data[,"outliers"]);
  }
  verbose && exit(verbose);

  len <- sapply(data, FUN=length);
  keep <- (len == nrow(map));
  data <- data[keep];
  data <- as.data.frame(data);

  data <- cbind(map, data);

  verbose && exit(verbose);

  data;
}, private=TRUE)



setMethodS3("updateDataFlat", "ChipEffectFile", function(this, data, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'data':
  names <- colnames(data);
  namesStr <- paste(names, collapse=", ");
  if (!"cell" %in% names)
    throw("Argument 'data' must contain a column 'cell': ", namesStr);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose2 <- -as.integer(verbose)-2;

  verbose && enter(verbose, "Storing flat data to file");

  # Encode?
  if ("outliers" %in% names) {
    data[,"outliers"] <- -as.integer(data[,"outliers"]);
  }

  names <- gsub("theta", "intensities", names);
  names <- gsub("sdTheta", "stdvs", names);
  names <- gsub("outliers", "pixels", names);
  colnames(data) <- names;

  verbose && enter(verbose, "Updating file");
  indices <- data[,"cell"];
  keep <- (names %in% c("intensities", "stdvs", "pixels"));
  data <- data[,keep];
  pathname <- getPathname(this);
  pathname <- Arguments$getWritablePathname(pathname);
  .updateCel(pathname, indices=indices, data, verbose=verbose2);
  verbose && exit(verbose);

  verbose && exit(verbose);
  invisible(data);
}, private=TRUE)



setMethodS3("mergeGroups", "ChipEffectFile", function(this, fcn, fields=c("theta", "sdTheta"), ..., pathname, overwrite=FALSE, verbose=FALSE) {
  # Argument 'fcn':
  if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", class(fcn)[1]);
  }

  # Argument 'pathname':
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=!overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Merging groups");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Test 'fcn':
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Testing merge function");
  for (size in 1:10) {
    y <- matrix(1000+1:(size*4), nrow=size);
    yOut <- fcn(y);
    if (!identical(dim(yOut), dim(y))) {
      throw("Function 'fcn' must not change the dimension of the data: ",
                                  paste(dim(yOut), collapse="x"), " != ",
                                            paste(dim(y), collapse="x"));
    }
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get flat (unit, group, cell) map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  map <- getUnitGroupCellMap(this, verbose=less(verbose));
  verbose && str(verbose, map);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge data for each unit size separately (in reverse order!)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get all possible unit sizes (number of groups per unit)
  uSizes <- sort(unique(data[,"group"]));
  verbose && cat(verbose, "Different number of groups per unit identified:");
  verbose && print(verbose, uSizes);

  data <- getDataFlat(this, units=map, ..., verbose=less(verbose));
  verbose && str(verbose, data);

  for (size in rev(uSizes)) {
    verbose && enter(verbose, "Unit size ", size);
    # Identify the units of that size
    idxs <- which(data[,"group"] == size);
    unitsS <- data[idxs, "unit"];

    # Get the subset of the data for such units
    idxs <- which(data[,"unit"] %in% unitsS);

    for (field in fields) {
      # Extract signals as a matrix where each column is one unit
      y <- data[idxs, field];
      y <- matrix(y, nrow=size);
      verbose && str(verbose, y);
      y <- fcn(y);
      verbose && str(verbose, y);

      # Update data table
      data[idxs, field] <- as.vector(y);
    }

    # Get the cells where to store the merged data
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Copy CEL file and update the copy
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Storing merged data");
  verbose && cat(verbose, "Pathname: ", pathname);

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing");
  cfN <- createFrom(this, filename=pathname, path=NULL, methods="create", clear=TRUE, verbose=less(verbose));
  verbose && print(verbose, cfN);
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing merged data");
  updateDataFlat(cfN, data=data, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && exit(verbose);


  verbose && exit(verbose);

  cfN;
}, protected=TRUE)


setMethodS3("extractMatrix", "ChipEffectFile", function(this, ..., field=c("theta", "sdTheta")) {
  # Argument 'field':
  field <- match.arg(field);

  NextMethod("extractMatrix", field=field);
})



############################################################################
# HISTORY:
# 2012-11-28
# o MEMORY: readUnits() for ChipEffectFile no longer caches results
#   by default.
# 2012-10-14
# o CLEANUP: createParamCdf() for ChipEffectFile no longer support
#   '<chipType>-monocell' filenames.  If detected, an informative
#   error is thrown.
# 2009-05-19
# o Now testing for file permissions for creat-/writ-/updating files/dirs.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: getUnitGroupCellMap().
# 2008-05-08
# o BUG FIX: getUnitGroupCellMap() of ChipEffectFile gave an error if
#   argument 'units' had zero length (non-NULL).
# 2008-04-21
# o getCellMap() is now defunct.
# 2008-03-11
# o Now getUnitGroupCellMap() of ChipEffectFile file caches smaller objects.
# 2008-02-28
# o Now a (unit,group,cell) map has class UnitGroupCellMap and no longer
#   ChipEffectFileCellMap.
# 2008-02-22
# o Added extractMatrix().
# o Renamed getCellMap() to getUnitGroupCellMap().
# o Now getCellMap() pass on the 'force' argument all the way.
# 2008-02-20
# o Now updateDataFile() only encodes the "outliers" field, if it part of
#   the input.  If the input is in "raw data", i.e. "pixels", it won't be
#   encoded.
# o Now fromDataFile() of ChipEffectFile no longer replicates the
#   "chipEffects" if already part of the filename.
# 2007-12-11
# o BUG FIX: getCellMap() of ChipEffectFile was broken.
# 2007-12-10
# o Now fromDataFile() of ChipEffectSet accepts argument 'cdf'.
# 2007-11-20
# o MEMORY OPTIMIZATION: Now getCellMap() builds data in chunks if
#   argument 'units' is NULL.
# 2007-09-12
# o Now getCellMap() of ChipEffectFile caches (large) results to file.
# 2007-08-16
# o Made findUnitsTodo() on ChipEffectFile much more memory efficient.
#   Before it could consume 1-2GB for the GenomeWideSNP_6 chip, but now
#   it is consuming ~100MB.  This was done by using lapplyInChunks().
# 2007-08-09
# o ChipEffectFile$fromDataFile() now creates CEL files with upper-case
#   filename extension "*.CEL", not "*.cel".  The reason for this is that
#   some software don't recognize lower case filename extensions :(
#   Note: The above modification is safe because the above method first
#   renames filename extensions in lower case to be in upper case, if
#   they exist.
# 2007-08-02
# o BUG FIX: getCellMap() would give 'Error in verbose && cat("Unique number
#   of groups per unit: ", paste(uUnitSizes,...', if verbose was on.
# 2007-07-19
# o Added more verbose output to getCellMap().
# 2007-06-11
# o BUG FIX: Called non-existing 'cf' instead of 'this' in mergeGroups()
#   of ChipEffectFile.  This function was never used anyway.
# 2007-02-20
# o Added mergeGroups().
# 2007-02-19 /HB
# o BUG FIX: getCellMap() did not handle units with different number of
#   groups.
# 2007-02-15 /KS
# o BUG FIX: getCellMap() did not handle units with other than one group.
# 2007-02-09
# o Updated the file cache sub directory.
# 2007-01-10
# o Now fromDataFile() looks for chip effect files named using the "sample
#   name" only (not tags) file format, and renames it to the full name
#   format.  This is a "patch" so we don't have to reestimate old data sets.
# 2007-01-09
# o Now fromDataFile() generates a file with the full name (name + tags) of
#   the input file and not just the name.
# 2007-01-07
# o TO DO: Speed up getCellMap(), e.g. using more clever file caching.
# 2007-01-06
# o findUnitsTodo() cache cell-index vector to file if all units are to
#   be scanned.
# 2007-01-05
# o Removed getSampleNames().
# 2007-10-02
# o TO DO: Static fromDataFile() does not really need argument 'df'.
# 2006-12-02
# o BUG FIX: getCellMap(..., units=NULL) did not work.
# 2006-11-28
# o Added trial version of updateDataFlat(). Seems to work. Will speed
#   up a few things.
# o Added trial versions of getCellMap() and getDataFlat().
# 2006-11-28
# o Added argument 'cache' to readUnits() to specify if the result should
#   be cached or not.
# o BUG FIX: The caching mechanism of readUnits() of ChipEffectFile was
#   not sensitive to the class of the ChipEffectFile object.  Added a
#   'class' element to the cache key.  Then each subclass must override
#   this method too.
# 2006-10-06
# o Now chip effect files use filename tag 'chipEffects' with a comma in
#   front (instead of suffix "-chipEffects').  This way getName() will
#   return the sample name without any extra endings.
# 2006-09-11
# o Great! Using the specially designed CDFs and CEL files to store
#   estimates is much faster and smaller than using the originally
#   structured CDF and CEL files.  Now storing the estimates takes a much
#   smaller part of the overall fitting algorithm.
# 2006-09-10
# o Starting to make use of specially design CDFs and CEL files for storing
#   chip effects.  This make getFirstCellIndices() obsolete.
# o Added createParamCdf().
# 2006-08-26
# o Created.  Have to store chip-effect estimates too.  Currently we use
#   the existing CEL/CDF structure for this, but those are unnecessarily
#   large for this.  Later this will be done in special CEL files with a
#   custom CDF file (possible virtual).  This requires that affxparser can
#   create empty CEL files from scratch, which is on the to-do list.
############################################################################
