###########################################################################/**
# @RdocClass AffymetrixCelFile
#
# @title "The AffymetrixCelFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixCelFile object represents a single Affymetrix CEL file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "aroma.core::AromaMicroarrayDataFile".}
#   \item{cdf}{An optional @see "AffymetrixCdfFile" making it possible to
#     override the default CDF file as specified by the CEL file header.
#     The requirement is that its number of cells must match that of
#     the CEL file.
#     If @NULL, the CDF structure is inferred from the the chip type
#     as specified in the CEL file header.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{For developers}{
#   If you subclass this class, please make sure to query the
#   @see "AffymetrixCdfFile" object (see @seemethod "getCdf") whenever
#   querying CDF information.  Do not use the CDF file inferred from the
#   chip type in CEL header, unless you really want it to be hardwired that
#   way, otherwise you will break to possibility to override the CDF
#   structure.
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically part of an @see "AffymetrixCelSet".
# }
#*/###########################################################################
setConstructorS3("AffymetrixCelFile", function(..., cdf=NULL) {
  this <- extend(AffymetrixFile(...), "AffymetrixCelFile",
    "cached:.header" = NULL,
    "cached:.lastPlotData" = NULL,
    .cdf = NULL
  )

  if (!is.null(cdf))
    setCdf(this, cdf);

  pathname <- getPathname(this);
  if (!is.null(pathname)) {
    # Make sure the name is non-empty
    name <- getName(this);
    if (!nzchar(name)) {
      throw("An ", class(this)[1], " must have a name of at least length one: ", pathname);
    }
  }

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("clone", "AffymetrixCelFile", function(this, ..., verbose=TRUE) {
  # Clone itself (and clear the cached fields)
  object <- NextMethod("clone", clear=TRUE);

  # Clone the CDF here.
  if (!is.null(object$.cdf))
    object$.cdf <- clone(object$.cdf);

  object;
}, protected=TRUE)


setMethodS3("getExtensionPattern", "AffymetrixCelFile", function(static, ...) {
  "[.](cel|CEL)$";
}, static=TRUE, protected=TRUE)


setMethodS3("getFileFormat", "AffymetrixCelFile", function(this, asString=TRUE, ...) {
  # Default
  ver <- NA_integer_
  verStr <- NA_character_

  pathname <- getPathname(this)
  if (isFile(pathname)) {
    # Read CEL header
    raw <- readBin(pathname, what=raw(), n=10)

    if (raw[1] == 59) {
      ver <- 1L
      verStr <- "v1 (binary; CC)"
    } else if (raw[1] == 64) {
      ver <- 4L
      verStr <- "v4 (binary; XDA)"
    } else if (rawToChar(raw[1:5]) == "[CEL]") {
      ver <- 3L
      verStr <- "v3 (text; ASCII)"
    }
  }

  if (asString) ver <- verStr

  ver
})

setMethodS3("as.character", "AffymetrixCelFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("File format: %s", getFileFormat(this, asString=TRUE)));
  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  s <- c(s, sprintf("Timestamp: %s", as.character(getTimestamp(this))));
  s;
}, protected=TRUE)


setMethodS3("getIdentifier", "AffymetrixCelFile", function(this, ..., force=FALSE) {
  identifier <- this$.identifier
  if (force || is.null(identifier)) {
    if (!isFile(this)) return(NA_character_)

    # Get header
    hdr <- getHeader(this);

    # AD HOC fix for backward compatibility.
    # Ideally we should not create an identifier based
    # on the pathname. /HB 2011-02-24
    pathname <- hdr$filename;
    path <- dirname(pathname);
    path <- dropRootPathTags(path, depth=2);
    filename <- basename(pathname);
    pathname <- file.path(path, filename);
    hdr$filename <- pathname;

    # Get subset of data
    nbrOfCells <- hdr$total;
    mid <- nbrOfCells %/% 2;
    subset <- seq(from=mid - 500, to=mid + 500);
    data <- getData(this, indices=subset);
    key <- list(hdr=hdr, data=data);
    id <- getChecksum(key);
    this$.identifier <- id;
  }
  identifier;
}, private=TRUE)




###########################################################################/**
# @RdocMethod fromFile
#
# @title "Defines an AffymetrixCelFile object from a CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of to the file.}
#  \item{path}{The path to the file.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelFile" object.
#  If the file is not found or if it is of the wrong file format, an
#  error is thrown.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromFile", "AffymetrixCelFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }


  # WORKAROUND: Currently the affxparser code crash R if the file is not
  # a valid CEL file.  The best we can do now is to test against the
  # filename.
  isCel <- (regexpr("[.](c|C)(e|E)(l|L)$", pathname) != -1);
  if (!isCel) {
    throw("Could not read CEL file. Filename format error: ", pathname);
  }

  # Create a new instance of the same class
  newInstance(static, pathname);
}, static=TRUE, protected=TRUE)





###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "AffymetrixCelFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    if (!isFile(this)) return(AffymetrixCdfFile())

    hdr <- getHeader(this);
    chipType <- hdr$chiptype;
    nbrOfCells <- hdr$total;
    cdf <- AffymetrixCdfFile$byChipType(chipType, nbrOfCells=nbrOfCells);
    this$.cdf <- cdf;
  }
  cdf;
})


setMethodS3("getUnitNamesFile", "AffymetrixCelFile", function(this, ...) {
  getCdf(this, ...);
})

setMethodS3("getUnitTypesFile", "AffymetrixCelFile", function(this, ...) {
  getCdf(this, ...);
})


###########################################################################/**
# @RdocMethod setCdf
#
# @title "Sets the CDF structure for this CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{...}{Not used.}
#   \item{.checkArgs}{(Internal) If @FALSE, arguments are not validated.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "getCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("setCdf", "AffymetrixCelFile", function(this, cdf, ..., .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'cdf':
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

    # Assure that the CDF is compatible with the CEL file
    if (nbrOfCells(cdf) != nbrOfCells(this)) {
      throw("Cannot set CDF. The specified CDF structure ('", getChipType(cdf), "') is not compatible with the chip type ('", getChipType(this), "') of the CEL file. The number of cells do not match: ", nbrOfCells(cdf), " != ", nbrOfCells(this));
    }

    # Nothing to do?
#    oldCdf <- getCdf(this);
#    if (equals(cdf, oldCdf))
#      return(invisible(this));
  }

  # Have to clear the cache
  clearCache(this);

  this$.cdf <- cdf;

  invisible(this);
})



###########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets the header of the CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @list structure as returned by @see "affxparser::readCelHeader".
#  If file does not exists, then @NULL is returned.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getHeader", "AffymetrixCelFile", function(this, ...) {
  header <- this$.header;
  if (is.null(header) && isFile(this)) {
    pathname <- getPathname(this);
    header <- .readCelHeader(pathname);
    this$.header <- header;
  }
  header;
}, private=TRUE)


setMethodS3("getHeaderV3", "AffymetrixCelFile", function(this, ...) {
  if (!isFile(this)) return(NULL)

  # Get the CEL header
  header <- getHeader(this);

  # Get the CEL v3 header
  header <- header$header;

  # Parse it
  header <- unlist(strsplit(header, split="\n"));
  header <- trim(header);
  header <- header[nzchar(header)];
  header <- strsplit(header, split="=");
  names <- sapply(header, FUN=function(s) s[1]);
  header <- lapply(header, FUN=function(s) s[-1]);
  names(header) <- names;

  header;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getTimestamp
#
# @title "Gets the timestamp in the CEL header"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{format}{The default format string for parsing the time stamp.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a \code{POSIXct} object.
#  The parsed string containing the timestamp is returned as
#  attribute \code{text}.
# }
#
# @author "HB"
#
# \seealso{
#   Internally, @see "base::strptime" is used to parse the time stamp.
#   @see "base::DateTimeClasses" for more information on \code{POSIXct}.
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getTimestamp", "AffymetrixCelFile", function(this, format="%m/%d/%y %H:%M:%S", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getTimestampFromDatHeader <- function(headerDAT, chipType) {
    # Find the element with a date. It is part of the same string as the
    # one containing the chip type.  Get the chip type from the header.
    pattern <- sprintf(" %s.1sq ", chipType);
    header <- grep(pattern, headerDAT, value=TRUE);

    # Fallback, in case a early-access chiptype string (or none) is used
    # in the DAT header.
    if (length(header) == 0L) {
      pattern <- " (.*).1sq ";
      header <- grep(pattern, headerDAT, value=TRUE);
    }

    hasTimestamp <- FALSE;
    if (length(header) == 1L) {
      # Extract the date timestamp
      pattern <- ".*([01][0-9]/[0-3][0-9]/[0-9][0-9] [0-2][0-9]:[0-5][0-9]:[0-5][0-9]).*";

      hasTimestamp <- (regexpr(pattern, header) != -1L);
    }

    if (hasTimestamp) {
      timestamp <- gsub(pattern, "\\1", header);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Alternative:
      # Could use a pattern, but if a different timestamp than the American is
      # used, this wont work.  Instead assume a fixed location.
      # From the DAT header specification (Affymetrix Data File Formats, April
      # 2006), we know that the date and the timestamp is 18 characters long.
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ##   nTemp <- 7;
    ##   nPower <- 4;
    ##   nTimestamp <- 18;
    ##   # Expected start position
    ##   pos <- nTemp + 1 + nPower + 1;
    ##   # ...however, the files we have start at the next position. /HB 2006-12-01
    ##   pos <- pos + 1;
    ##   timestamp <- substring(header, first=pos, last=pos+nTimestamp-1);

      timestamp <- trim(timestamp); # Unnecessary?

      # Parse the identified timestamp into POSIXct
      res <- strptime(timestamp, format=format, ...);
    } else {
      res <- NULL;
    }

    # If no valid timestamp was found, return NA.
    if (length(as.character(res)) == 0) {
      res <- as.POSIXct(NA);
    }

    # Keep the non-parsed timetstamp etc for debugging.
    attr(res, "text") <- timestamp;
    attr(res, "header") <- header;
    attr(res, "headerDAT") <- headerDAT;

    res;
  } # getTimestampFromDatHeader()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'format':
  format <- Arguments$getCharacter(format, length=c(1,1));

  if (!isFile(this)) return(as.POSIXct(NA))

  chipType <- getHeader(this)$chiptype;

  fileFormat <- getFileFormat(this, asString=FALSE);

  if (fileFormat == 1) {
    suppressWarnings({
      hdr <- .readCcgHeader(getPathname(this));
    });

    # Get the DAT header
    params <- hdr$dataHeader$parents[[1]]$parameters;
    headerDAT <- params[["affymetrix-dat-header"]];
    if (is.null(headerDAT)) {
      headerDAT <- params[["affymetrix-partial-dat-header"]];
    }

    # Extract the timestamp
    res <- getTimestampFromDatHeader(headerDAT, chipType=chipType);
  } else if (fileFormat %in% c(3,4)) {
    # Get the DAT header
    header <- getHeaderV3(this);
    headerDAT <- header$DatHeader;

    # Extract the timestamp
    res <- getTimestampFromDatHeader(headerDAT, chipType=chipType);
  } else {
    res <- as.POSIXct(NA);
  }

  res;
}, private=TRUE)



setMethodS3("nbrOfCells", "AffymetrixCelFile", function(this, ...) {
  hdr <- getHeader(this)
  if (!is.list(hdr)) return(NA_integer_)
  hdr$total
})


###########################################################################/**
# @RdocMethod getChipType
#
# @title "Gets the chip type for this CEL file"
#
# \description{
#  @get "title" \emph{according} to the @see "AffymetrixCdfFile" object.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getChipType", "AffymetrixCelFile", function(this, ...) {
  unf <- getUnitNamesFile(this);
  if (!isFile(unf)) return(NA_character_)
  getChipType(unf, ...);
}, private=TRUE)


setMethodS3("getPlatform", "AffymetrixCelFile", function(this, ...) {
  "Affymetrix";
}, private=TRUE)



###########################################################################/**
# @RdocMethod readUnits
#
# @title "Reads CEL data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probeset).
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{cdf}{An alternative CDF structure to be used.  This overrides
#     the \code{units} arguments.}
#   \item{...}{Arguments passed to \code{getUnits()} of the
#     @see "AffymetrixCdfFile" class (if \code{cdf} was not specified),
#     but also to the @see "affxparser::readCelUnits" methods.}
# }
#
# \value{
#  Returns the @list structure that @see "affxparser::readCelUnits" returns.
# }
#
# \section{Caching}{
#   CEL data is neither cached in memory nor on file by this method.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "updateUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "AffymetrixCelFile", function(this, units=NULL, cdf=NULL, ..., stratifyBy=NULL, force=FALSE, cache=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve CDF structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(cdf)) {
    suppressWarnings({
      cdf <- readUnits(getCdf(this), units=units, stratifyBy=stratifyBy);
    });
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- getPathname(this);
  suppressWarnings({
    res <- .readCelUnits(pathname, cdf=cdf, dropArrayDim=TRUE, ...);
  })

  res;
}, private=TRUE)


###########################################################################/**
# @RdocMethod updateUnits
#
# @title "Updates CEL data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units.
# }
#
# @synopsis
#
# \arguments{
#   \item{data}{A @list structure consisting of grouped data similar to
#      what @seemethod "readUnits" returns.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns the @list structure that @see "affxparser::updateCelUnits" returns.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "updateUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("updateUnits", "AffymetrixCelFile", function(this, data, ...) {
  pathname <- getPathname(this);
  pathname <- Arguments$getWritablePathname(pathname);
  .updateCelUnits(pathname, data=data, ...);
}, private=TRUE)



###########################################################################/**
# @RdocMethod clearData
#
# @title "Clears all or a subset of the fields in a CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fields}{A @character @vector of fields to be cleared.}
#   \item{value}{A @numeric value to be written over the data.}
#   \item{...}{Not used.}
#   \item{.forSure}{If not @TRUE, an exception is thrown asking if the
#      method was called by mistake.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns (invisibly) the names of the fields cleared.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
#   Internally, @see "affxparser::updateCel" is used.
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("clearData", "AffymetrixCelFile", function(this, fields=c("intensities", "stdvs", "pixels"), value=0, ..., .forSure=FALSE, verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fields':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  if (!identical(.forSure, TRUE))
    throw("Did you call clearData() by mistake? If not, use .forSure=TRUE.");

  # Nothing do to?
  if (length(fields) == 0L) {
    verbose && cat(verbose, "No fields to be cleared.");
    return(invisible(fields));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clear
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Clearing Affymetrix CEL file");

  # Asserting file permissions
  pathname <- getPathname(this);
  pathname <- Arguments$getWritablePathname(pathname);

  verbose && cat(verbose, "Fields to be cleared: ",
                                             paste(fields, collapse=", "));
  bfr <- rep(value, length.out=nbrOfCells(this));
  intensities <- stdvs <- pixels <- NULL;
  if ("intensities" %in% fields)
    intensities <- bfr;
  if ("stdvs" %in% fields)
    stdvs <- bfr;
  if ("pixels" %in% fields)
    pixels <- bfr;

  .updateCel(pathname, intensities=bfr, stdvs=bfr, pixels=bfr);
  verbose && exit(verbose);

  invisible(fields);
}, static=TRUE, private=TRUE)




###########################################################################/**
# @RdocMethod readRawData
# @aliasmethod getData
#
# @title "Gets all or a subset of the fields in a CEL file"
#
# \description{
#  @get "title" for all or a subset of the cells.
# }
#
# @synopsis
#
# \arguments{
#   \item{indices}{A @numeric @vector of cell indices.  If @NULL, all cells
#     are considered.}
#   \item{fields}{Names of fields to be retrieved.}
#   \item{...}{Additional arguments passed to @see "affxparser::readCel".}
#   \item{drop}{If @TRUE and a single field is returned, then data is
#     returned as a @vector, otherwise as a @data.frame.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @data.frame of the fields requested (unless dimension dropped).
# }
#
# \section{Caching}{
#   Neither in-memory nor on-file caching is done by this method.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readRawData", "AffymetrixCelFile", function(this, indices=NULL, fields=c("xy", "intensities", "stdvs", "pixels"), ..., drop=FALSE, verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  readCel <- affxparser::readCel


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indices':
  nbrOfCells <- nbrOfCells(getCdf(this));
  if (is.null(indices)) {
  } else {
    indices <- Arguments$getIndices(indices, max=nbrOfCells, disallow="NaN");
    nbrOfCells <- length(indices);
  }

  # Argument 'fields':
  fields <- match.arg(fields, several.ok=TRUE);
  if (length(fields) == 0) {
    throw("Argument 'fields' is empty.");
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Special case: requesting zero indices?
  readZeroElements <- (length(indices) == 0 && !is.null(indices));
  if (readZeroElements) {
    # A work around...
    indices <- 1;
  }

  # Workaround for readCel() not handling NA indices
  if (!is.null(indices)) {
    nas <- which(is.na(indices));
    hasNAs <- length(nas);
    if (hasNAs)
      indices[nas] <- 1;
  } else {
    hasNAs <- FALSE;
  }

  cVerbose <- -(as.numeric(verbose) + 50);
  pathname <- getPathname(this);
  args <- list(
    filename=pathname,
    indices=indices,
    readHeader=FALSE,
    readIntensities=is.element("intensities", fields),
    readStdvs=is.element("stdvs", fields),
    readPixels=is.element("pixels", fields),
    readXY=is.element("xy", fields),
    readOutliers=FALSE,
    readMasked=FALSE,
    ...,
    verbose=cVerbose
  );

  fcn <- get("readCel", mode="function");
  keep <- intersect(names(args), names(formals(fcn)));
  args <- args[keep];
  cel <- do.call(readCel, args=args);

  # Sanity check
  stopifnot(is.list(cel));
  stopifnot(length(cel) > 0);

  if (hasNAs) {
    for (kk in seq_along(cel)) {
      naValue <- NA;
      storage.mode(naValue) <- storage.mode(cel[[kk]]);
      cel[[kk]][nas] <- naValue;
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clean up
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Split (x,y)
  isXY <- which(fields == "xy");
  if (length(isXY) > 0) {
    fields <- as.list(fields);
    fields[[isXY]] <- c("x", "y");
    fields <- unlist(fields, use.names=TRUE);
  }

  # Keep only requested fields
  if (!identical(names(cel), fields)) {
    cel <- cel[fields];

    # Sanity check
    stopifnot(is.list(cel));
    stopifnot(length(cel) > 0);
  }

  if (readZeroElements) {
    cel <- lapply(cel, FUN=.subset, integer(0));
  }

  # Drop dimensions?
  if (drop && length(cel) == 1) {
    cel <- cel[[1]];
  } else {
    # Return as data frame
    attr(cel, "row.names") <- seq_len(length(cel[[1]]));
    class(cel) <- "data.frame";
  }

  cel;
}, private=TRUE)


setMethodS3("range", "AffymetrixCelFile", function(this, ..., na.rm=TRUE) {
  x <- getData(this, ...);
  range(x, na.rm=na.rm);
}, protected=TRUE)


setMethodS3("readRawDataRectangle", "AffymetrixCelFile", function(this, xrange=c(0,Inf), yrange=c(0,Inf), fields=c("intensities", "stdvs", "pixels"), ..., drop=FALSE) {
  pathname <- getPathname(this);
  data <- .readCelRectangle(pathname, xrange=xrange, yrange=yrange, readIntensities=("intensities" %in% fields), readStdvs=("stdvs" %in% fields), readPixels=("pixels" %in% fields), readHeader=FALSE, readOutliers=FALSE, readMasked=FALSE);

  if (drop && length(data) == 1) {
    data <- data[[1]];
  }

  data;
}, private=TRUE)


setMethodS3("getRectangle", "AffymetrixCelFile", function(this, ...) {
  readRawDataRectangle(this, ...);
}, private=TRUE)



############################################################################
# HISTORY:
# 2012-12-10
# o BUG FIX: getTimestamp() for AffymetrixCelFile would throw an error
#   if the CEL file header did not have a timestamp.
# 2012-11-20
# o CLEANUP: Deprecated "[" and "[[", because they should be used to
#   subset files and not units.
# 2011-11-18
# o ROBUSTNESS: Added validiation of argument 'fields' to readRawData()
#   of AffymetrixCelFile and more internal sanity checks in that method.
# 2011-08-31
# o ROBUSTNESS: getTimestamp() for AffymetrixCelFile would throw "Error
#   in if (hasTimestamp) { : argument is of length zero" if the CEL file
#   had a DAT header with a non-standard chip type string, e.g. an
#   early-access label or no label at all.  Updated the local/inner
#   getTimestampFromDatHeader() of getTimestamp() to also handle such
#   CEL files.  Thanks Irina Ostrovnaya at MSKCC for reporting on this.
# 2011-02-24
# o BACKWARD COMPATILITY: getIdentifier() for AffymetrixCelFile generates
#   a checksum id based on the CEL file header among other things.  Part
#   of this header information is the relative pathname of the CEL file,
#   which means that the identifier will be different depending on which
#   directory it lives in.  Ideally we'll exclude the pathname from this.
#   However, in the meanwhile we simply drop any tags from the root path
#   such that it is compatible with earlier version of aroma.*.
# 2011-02-22
# o ROBUSTNESS: getTimestamp() for AffymetrixCelFile no longer assumes
#   that the file's DAT header contains a timestamp and tries to translate.
#   Instead it first tests for the timestamp pattern, and returns NA
#   if not found.
# 2009-07-08
# o Added getUnitTypesFile() for AffymetrixCelFile.
# 2009-05-19
# o Now testing for file permissions before trying to update a CEL file.
# 2008-06-25
# o BUG FIX: getChipType() of AffymetrixCelFile did not pass down '...'
#   causing for instance getChipType(..., fullname=FALSE) to still return
#   tags for ChipEffectFile:s.
# 2008-05-09
# o Added getUnitNamesFile().
# o Added getPlatform().
# o Now AffymetrixCelFile inherits from AromaMicroarrayDataFile.
# 2008-05-08
# o BUG FIX: readRawData() did not handle a zero-length 'indices' argument;
#   it was interpreted as NULL, i.e. read everything.
# 2008-03-14
# o Renamed getRectangle() to readRawDataRectangle().
# 2008-03-13
# o Now readRawData(), formely known as getData(), handle NA indices and
#   drops the dimension if only one field and argument drop=TRUE.
# 2008-03-05
# o Now setCdf() also reports on the two incompatible chip types involved
#   if trying to set a CDF that is not compatible with a CEL file.
# 2008-02-21
# o Updated getTimestamp() to also support Calvin files containing only
#   a partial DAT header.
# 2008-01-30
# o Now getTimestamp() for AffymetrixCelFile works for Calvin files too.
# o Now getFileFormat() for AffymetrixCelFile takes argument 'asString'.
# 2007-07-09
# o Added getFileFormat() to AffymetrixCelFile.  This is also reported
#   by the print() method.
# 2007-05-09
# o BUG FIX: If no valid timestamp was identified in the CEL header by
#   getTimestamp() of AffymetrixCelFile, then as.character() would give
#   "Error in sprintf(fmt, ...) : zero-length argument".  Now it returns
#   NA instead as a fall back if no valid timestamp is found.
# 2007-03-23
# o BUG FIX: Called non-existing readHeader() instead of getHeader().
# 2007-03-05
# o Added setAttributesByTags().
# o Added setAttributeXY(), getAttributeXY(), and hasAttributeXY().
# 2007-02-12
# o Now getData() is using do.call() because it is faster. Unused arguments
#   are still ignored.
# 2007-02-04
# o Now getData() is call readCel() using doCall() so that unused arguments
#   in '...' are ignored.
# 2007-02-03
# o BUG FIX: getTimestamp() assumed a fix location in the CEL v3 header,
#   but that did not work for dChip exported CEL files.  Now, a US date
#   pattern is assumed and searched for.
# 2007-01-12 [KS]
# o Moved image270() and writeSpatial() to AffymetrixCelFile.PLOT.R.
# 2006-12-18 [KS]
# o Add "takeLog" argument (logical) to image270.  If true, take the log2
#   before plotting.  Can be more informative than natural scale.
# 2006-12-14
# o Removed getSampleName() which gives the same as getName().
# 2006-12-11
# o Now the timestamp is also reported for singel CEL files.
# o BUG FIX: getHeaderV3() would throw an error if there was an empty V3
#   header fields.  This was the reason why getTimestamp() gave an error
#   on some 100K chips.
# 2006-12-01
# o Added getTimestamp().
# 2006-11-28
# o Arguments 'force' and 'cache' has to be in readUnits() to avoid being
#   passed from calls of subclasses.
# 2006-10-23
# o Update default value for argument 'fields' in getData().
# 2006-10-22
# o In order to speed up fromFile(), the CEL header is not read anymore.
# 2006-10-06
# o make sure cdf association is inherited
# 2006-08-28
# o Renamed getFields() to getData() because getFields() is "reserved"
#   for use in the Object class.
# 2006-08-27
# o Added nbrOfCells() because it is so common.
# o Added createFrom() which utilizes new functions copyFile() and
#   clearData(). It is also no longer static. This is more generic and
#   cleaner.  The new clearData() does also not require the CDF file
#   (in case that should be missing).
# 2006-08-25
# o Renamed getIntensities() to getFields() which returns a data frame.
# o Added image270() and writeSpatial().
# o Added methods "[" and "[[" mapped to readUnits().
# 2006-08-24
# o Added the option to specify an 'cdf' object, making it possible to
#   override the default CDF file according to the CEL header.  It is
#   important that all methods queries the AffymetrixCdfFile object
#   from getCdf() and not the one through the CEL header!
# o Added most Rdoc comments.
# 2006-07-21
# o Added readUnits().
# 2006-07-05
# o BUG FIX/WORKAROUND: Currently the affxparser code crash R if the file
#   is not a valid CEL file.  The best we can do now is to test that the
#   filename has suffix *.CEL.
# 2006-05-30
# o Added fromFile().
# 2006-03-30
# o Updated according to affxparser.
# 2006-03-23
# o Moved all SNP related methods into the new class AffymetrixSnpCelFile.
# 2006-03-18
# o Made probe indices one-based.
# 2006-03-04
# o Added support for remapping in readIntensities().  This is currently
#   not used for CEL files (only APD files), but was added for the future.
# 2006-03-02
# o Created.
############################################################################
