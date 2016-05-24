###########################################################################/**
# @RdocClass FirmaFile
#
# @title "The FirmaFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents scores calculated by the FIRMA algorithm.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixCelFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS, HB"
#
# \seealso{
#   An object of this class is typically part of a @see "FirmaSet".
# }
#
#*/###########################################################################
setConstructorS3("FirmaFile", function(...) {
  this <- extend(ParameterCelFile(...), "FirmaFile");

  setEncodeFunction(this, function(groupData, ...) {
    groupData;
  })

  setDecodeFunction(this, function(groupData, ...) {
    groupData;
  })

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("findUnitsTodo", "FirmaFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Identifying non-assigned units in FIRMA file");

  verbose && cat(verbose, "Pathname: ", getPathname(this));
  if (is.null(units)) {
    units <- 1:nbrOfUnits(getCdf(this));
  }

  verbose && exit(verbose);

  # Read pixels from each unit
  verbose && enter(verbose, "Reading data for these ", length(units), " units");
#  value <- .readCelUnits(getPathname(this), units=units, readIntensities=FALSE,
#                        readStdvs=FALSE, readPixels=TRUE, dropArrayDim=TRUE);

  value <- readUnits(this, units=units, readIntensities=FALSE,
                        readStdvs=FALSE, readPixels=TRUE, force=force);

  # Identify units for which all pixels == 0.

  allZeroPixels <- sapply(value, function(x) {all(x[[1]][[1]]==0)}, USE.NAMES=FALSE);

  value <- which(allZeroPixels);
  if (!is.null(units))
    value <- units[value];
  verbose && cat(verbose, "Looking for pixels == 0 indicating non-assigned units:");
  verbose && str(verbose, value);

  verbose && exit(verbose);

  value;
})


setMethodS3("createParamCdf", "FirmaFile", function(static, sourceCdf, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Creating CDF for FIRMA results");
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
    verbose && cat(verbose, "Will create CDF for the FIRMA results files from the
original CDF. NOTE: This will take several minutes or more!");
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



setMethodS3("readUnits", "FirmaFile", function(this, units=NULL, cdf=NULL,
..., force=FALSE, cache=FALSE, verbose=FALSE) {
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
    verbose && cat(verbose, "readUnits.FirmaFile(): Returning cached data");
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
    verbose && cat(verbose, "readUnits.FirmaFile(): Updating cache");
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
setMethodS3("getCellIndices", "FirmaFile", function(this, ..., .cache=TRUE) {
  cdf <- getCdf(this);
  getCellIndices(cdf, ...);
}, protected=TRUE)




setMethodS3("updateUnits", "FirmaFile", function(this, units=NULL, cdf=NULL, data, ...) {
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  NextMethod("updateUnits", cdf=cdf, data=data);
}, private=TRUE);


setMethodS3("fromDataFile", "FirmaFile", function(static, df=NULL, filename=sprintf("%s,FIRMAscores.CEL", getFullName(df)), path, name=getName(df), cdf=NULL, ..., verbose=FALSE) {
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



  # Rename lower-case *.cel to *.CEL, if that is the case.  Old versions
  # of the package generated lower-case CEL files. /HB 2007-08-09
  pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);

  if (!isFile(pathname)) {
    verbose && enter(verbose, "Creating FIRMA results file");
    verbose && cat(verbose, "Pathname: ", pathname);

    # Get CDF for chip effects
    if (is.null(cdf)) {
      cdf <- createParamCdf(static, getCdf(df), verbose=less(verbose));
    }

    # Get CDF header
    cdfHeader <- getHeader(cdf);

    # Build a valid CEL header
    celHeader <- .cdfHeaderToCelHeader(cdfHeader, sampleName=name);

    # Add some extra information about what the CEL file is for
    params <- c(Description="This CEL file contains FIRMA results calculated by the aroma.affymetrix package.");
    parameters <- gsub(" ", "_", params);
    names(parameters) <- names(params);
    parameters <- paste(names(parameters), parameters, sep=":");
    parameters <- paste(parameters, collapse=";");
    parameters <- paste(celHeader$parameters, parameters, "", sep=";");
    parameters <- gsub(";;", ";", parameters);
    parameters <- gsub(";$", "", parameters);
    celHeader$parameters <- parameters;

    # Write to a temporary file
    pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

    # Create the CEL file
    .createCel(pathnameT, header=celHeader, ..., verbose=less(verbose));

    # Rename temporary file
    popTemporaryFile(pathnameT, verbose=verbose);

    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Defining FIRMA results file");
  verbose && cat(verbose, "Pathname: ", pathname);
  res <- newInstance(static, pathname);
  # Inherit the CDF?
  if (!is.null(cdf))
    setCdf(res, cdf);
  verbose && exit(verbose);

  res;
}, static=TRUE, private=TRUE)


setMethodS3("getUnitGroupCellMap", "FirmaFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (inherits(units, "UnitGroupCellMap")) {
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
# BUG!  Fix this in ChipEffectFile.R
#  unitSizes <- unlist(lapply(cells, FUN=length), use.names=FALSE);
  unitSizes <- unlist(lapply(cells, FUN=function(unit){
    length(.subset2(unit,"groups"));
  }), use.names=FALSE);
  cells <- unlist(cells, use.names=FALSE);
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating return data frame");
  uUnitSizes <- unique(unitSizes);
  if (is.null(units)) {
    cdf <- getCdf(this);
    units <- seq_len(nbrOfUnits(cdf));
  }

# BUG!  Fix this in ChipEffectFile.R
#  units <- rep(units, each=unitSizes);
  units <- rep(units, unitSizes);

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

  class(map) <- c("UnitGroupCellMap", class(map));

  map;
}, private=TRUE)


setMethodS3("getDataFlat", "FirmaFile", function(this, units=NULL, fields=c("intensities", "stdvs", "pixels"), ..., verbose=FALSE) {
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
  celFields <- c(intensities="intensities", stdvs="stdvs", pixels="pixels");
  suppressWarnings({
    data <- getData(this, indices=map[,"cell"], fields=celFields[fields]);
  })
  rownames(data) <- seq_len(nrow(data));  # Work around?!? /HB 2006-11-28

  verbose && str(verbose, data);

  verbose && exit(verbose);

  len <- sapply(data, FUN=length);
  keep <- (len == nrow(map));
  data <- data[keep];
  data <- as.data.frame(data);
  data <- cbind(map, data);

  verbose && exit(verbose);

  data;
}, private=TRUE)


setMethodS3("updateDataFlat", "FirmaFile", function(this, data, ..., verbose=FALSE) {
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

  colnames(data) <- names;

  verbose && enter(verbose, "Updating file");
  indices <- data[,"cell"];
  keep <- (names %in% c("intensities", "stdvs", "pixels"));
  data <- data[,keep];
  pathname <- getPathname(this);
  pathname <- Arguments$getWritablePathname(pathname):
  .updateCel(pathname, indices=indices, data, verbose=verbose2);
  verbose && exit(verbose);

  verbose && exit(verbose);
  invisible(data);
}, private=TRUE)


setMethodS3("extractMatrix", "FirmaFile", function (this, ..., field=c("intensities", "stdvs", "pixels")) {
  # Argument 'field':
  field <- match.arg(field);

  NextMethod("extractMatrix", field=field);
})


############################################################################
# HISTORY:
# 2012-10-14
# o CLEANUP: createParamCdf() for FirmaFile no longer support
#   '<chipType>-monocell' filenames.  If detected, an informative
#   error is thrown.
# 2010-05-12
# o ROBUSTNESS: When fromDataFile() of FirmaFile creates a file, it
#   is created first as a temporary file which is then renamed.  This
#   lowers the risk of generating incomplete chip-effect files.
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
# 2008-02-28 [HB]
# o Now a (unit,group,cell) map has class UnitGroupCellMap and no longer
#   ChipEffectFileCellMap.
# 2008-02-22 [HB]
# o Added extractMatrix().
# 2007-08-09
# o FirmaFile$fromDataFile() now creates CEL files with upper-case
#   filename extension "*.CEL", not "*.cel".  The reason for this is that
#   some software don't recognize lower case filename extensions :(
# 2007-02-09
# o Created (based on ChipEffectFile.R).
############################################################################
