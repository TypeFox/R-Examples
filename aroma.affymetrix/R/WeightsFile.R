###########################################################################/**
# @RdocClass WeightsFile
#
# @title "The WeightsFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents weights calculated from residuals of
#  probe-level models.
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
# @author "HB, KS"
#
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getWeightsSet()} method for the @see "ProbeLevelModel" class.
#   An object of this class is typically part of a @see "WeightsSet".
# }
#*/###########################################################################
setConstructorS3("WeightsFile", function(..., probeModel=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  this <- extend(ParameterCelFile(...), "WeightsFile",
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )

  setEncodeFunction(this, function(groupData, ...) {
    wts <- .subset2(groupData, "wts");
    wtsStdvs <- .subset2(groupData, "wtsStdvs");

    res <- list();
    if (!is.null(wts))
      res$intensities <- wts;
    if (!is.null(wtsStdvs))
      res$stdvs <- wtsStdvs;

    res;
  })

  setDecodeFunction(this, function(groupData, ...) {
    res <- list();
    if (!is.null(groupData$intensities))
      res$wts <- groupData$intensities;
    if (!is.null(groupData$stdvs))
      res$wtsStdvs <- groupData$stdvs;
    res;
  })

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("as.character", "WeightsFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  s;
}, protected=TRUE)


setMethodS3("getParameters", "WeightsFile", function(this, ...) {
  params <- NextMethod("getParameters");
  params$probeModel <- this$probeModel;
  params;
}, protected=TRUE)


setMethodS3("fromDataFile", "WeightsFile", function(static, df=NULL, filename=sprintf("%s,weights.CEL", getFullName(df)), path, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  if (!is.null(df)) {
    df <- Arguments$getInstanceOf(df, "AffymetrixCelFile");
  }

  # Argument 'cdf':
  if (is.null(cdf)) {
    if (is.null(df))
      throw("Either argument 'df' or 'cdf' must specified.");
  } else {
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");
  }

  # Argument 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Backward compatibility patch for now. Before residual files
  # only carried on the sample name, but not the tags. If such a
  # file is detected, it is renamed.
  # This should be removed in future versions. /HB 2007-01-10
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Pathname: ", pathname);
  res <- createFrom(df, filename=pathname, methods="create", clear=TRUE, ...);

  # Don't forget to return a ResidualFile object
  res <- fromFile(static, filename=pathname, verbose=less(verbose));
  # Inherit the CDF?
  if (!is.null(cdf))
    setCdf(res, cdf);
  verbose && print(verbose, res);

  res;
}, static=TRUE, private=TRUE)



setMethodS3("readUnits", "WeightsFile", function(this, units=NULL, cdf=NULL, ..., force=FALSE, cache=FALSE, verbose=FALSE) {
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
    verbose && cat(verbose, "readUnits.WeightsFile(): Returning cached data");
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
    verbose && cat(verbose, "readUnits.WeightsFile(): Updating cache");
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
setMethodS3("getCellIndices", "WeightsFile", function(this, ..., .cache=TRUE) {
  cdf <- getCdf(this);
  getCellIndices(cdf, ...);
}, protected=TRUE)


setMethodS3("updateUnits", "WeightsFile", function(this, units=NULL, cdf=NULL, data, ...) {
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  NextMethod("updateUnits", cdf=cdf, data=data);
}, private=TRUE);



setMethodS3("findUnitsTodo", "WeightsFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Identifying non-calculated units in residual file");

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

    verbose && enter(verbose, "Reading CDF cell indices");
    idxs <- getCellIndices(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);

    verbose && enter(verbose, "Extracting first cell in the first block for each unit");
    idxs <- .applyCdfGroups(idxs, function(groups) {
      # == groups[[1]]$indices[1];
      .subset(.subset2(.subset2(groups, 1), "indices"), 1)
    });
    verbose && exit(verbose);
    idxs <- unlist(idxs, use.names=FALSE);

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



setMethodS3("getUnitGroupCellMap", "WeightsFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (inherits(units, "WeightsFileCellMap")) {
    return(units);
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


  # Is 'units' already a CDF list?
  if (is.list(units)) {
    # No fancy validation for now.
    cells <- units;
    cdf <- getCdf(this);
    units <- indexOf(cdf, names=names(units));
    if (any(is.na(units))) {
      throw("Argument 'units' is of unknown structure.");
    }
    verbose && enter(verbose, "Argument 'cells' is already a CDF cell-index structure");
  } else {
    verbose && enter(verbose, "Retrieving cell indices for specified units");
    # Get the cells to read
    cells <- getCellIndices(this, units=units, force=force, verbose=less(verbose));
  }

  unitNames <- names(cells);
  unitSizes <- unlist(lapply(cells, FUN=length), use.names=FALSE);
  cells <- unlist(cells, use.names=FALSE);
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating return data frame");
  uUnitSizes <- unique(unitSizes);
  if (is.null(units)) {
    cdf <- getCdf(this);
    units <- seq_len(nbrOfUnits(cdf));
  }
  units <- rep(units, each=unitSizes);

  # The following is too slow:
  #  groups <- sapply(unitSizes, FUN=function(n) seq_len(n));

  # Instead, updated size by size
  naValue <- as.integer(NA);
  groups <- matrix(naValue, nrow=max(uUnitSizes), ncol=length(unitNames));
  for (size in uUnitSizes) {
    cc <- which(unitSizes == size);
    seq <- seq_len(size);
    groups[seq,cc] <- seq;
  }
  groups <- groups[!is.na(groups)];
  map <- data.frame(unit=units, group=groups, cell=cells);
  verbose && exit(verbose);

  verbose && exit(verbose);

  class(map) <- c("WeightsFileCellMap", class(map));

  map;
}, private=TRUE)



setMethodS3("getDataFlat", "WeightsFile", function(this, units=NULL, fields=c("theta", "sdTheta", "outliers"), ..., verbose=FALSE) {
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



setMethodS3("updateDataFlat", "WeightsFile", function(this, data, ..., verbose=FALSE) {
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

  # Encode
  names <- gsub("theta", "intensities", names);
  names <- gsub("sdTheta", "stdvs", names);
  names <- gsub("outliers", "pixels", names);
  colnames(data) <- names;
  if ("pixels" %in% names) {
    data[,"pixels"] <- -as.integer(data[,"pixels"]);
  }

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


setMethodS3("getImage", "WeightsFile", function(this, zrange=c(-1,1)*15, transform=log2, palette=rainbow(256), ...) {
  NextMethod("getImage", zrange=zrange, transform=transform, palette=palette);
})

setMethodS3("writeImage", "WeightsFile", function(this, ..., tags=c("*", "log2", "rainbow")) {
  NextMethod("writeImage", tags=tags);
})


############################################################################
# HISTORY:
# 2012-11-28
# o MEMORY: readUnits() for WeightsFile no longer cache results by default.
# 2009-05-19
# o Now testing for file permissions for creat-/writ-/updating files/dirs.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: getUnitGroupCellMap().
# 2008-05-08
# o BUG FIX: getUnitGroupCellMap() gave an error if argument 'units' had
#   zero length (non-NULL).
# 2008-04-21
# o getCellMap() is now defunct.
# 2007-08-09
# o WeightFile$fromDataFile() now creates CEL files with upper-case
#   filename extension "*.CEL", not "*.cel".  The reason for this is that
#   some software don't recognize lower case filename extensions :(
# 2007-04-12
# o BUG FIX: fromDataFile() of ResidualFile returned an AffymetrixCelFile
#   but not a ResidualFile.  This caused getResidualSet() of ProbeLevelModel
#   to return a ResidualSet containing AffymetrixCelFile:s.  The same
#   bug was found for the WeightFile class.  This problem was reported on
#   the mailing list on 2007-04-06.
# 2007-02-15
# o Created from ResidualFile.R.
############################################################################
