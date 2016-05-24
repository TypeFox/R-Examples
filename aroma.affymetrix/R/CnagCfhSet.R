###########################################################################/**
# @RdocClass CnagCfhSet
#
# @title "The CnagCfhSet class"
#
# \description{
#  @classhierarchy
#
#  An CnagCfhSet object represents a set of CNAG CFH files
#  with \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "CnagCfhFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \seealso{
#   @see "CnagCfhFile".
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("CnagCfhSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    lapply(files, FUN=function(df) {
      df <- Arguments$getInstanceOf(df, "CnagCfhFile", .name="files");
    })
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  extend(GenericDataFileSet(files=files, ...), "CnagCfhSet",
    "cached:.fileSize" = NULL
  )
})


setMethodS3("clone", "CnagCfhSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Cloning CNAG CFH set");

  # Clone itself and the files.  The call below will clear the cache!
  object <- NextMethod("clone", clear=TRUE, verbose=less(verbose));
  clearCache(object);

  # Clone the CDF (this will update the CDF of all file object)
  verbose && enter(verbose, "Cloning CDF");
  cdf <- clone(getCdf(object));
  verbose && exit(verbose);
  verbose && enter(verbose, "Adding CDF to CFH set");
  setCdf(object, cdf, .checkArgs=FALSE);
  verbose && exit(verbose);

  verbose && exit(verbose);

  object;
}, protected=TRUE)


setMethodS3("append", "CnagCfhSet", function(this, other, clone=TRUE, ..., verbose=FALSE) {
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, class(this)[1]);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Appending CFH set");
  verbose && print(verbose, other);

  # Validate chip type
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  for (file in getFiles(other)) {
    oCdf <- getCdf(file);
    oChipType <- getChipType(oCdf);
    if (!identical(oChipType, chipType)) {
      throw("Argument 'other' contains a CFH file of different chip type: ",
                                                oChipType, " != ", chipType);
    }
  }

  # Append other
  this <- NextMethod("append", other=other, clone=clone);

  # Set the same CDF for all CFH files
  verbose && enter(verbose, "Updating the CDF for all files");
  setCdf(this, cdf);
  verbose && exit(verbose);

  verbose && exit(verbose);

  this;
}, protected=TRUE)




###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the CNAG CFH set"
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
#  Returns a @character string.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "CnagCfhSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  s <- c(s, sprintf("Tags: %s", tags));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  n <- length(this);
  s <- c(s, sprintf("Number of arrays: %d", n));
  names <- getNames(this);
  s <- c(s, sprintf("Names: %s [%d]", hpaste(names), n));
  # Get CFH header timestamps
  ts <- getTimestamps(this);
  # Note: If ts <- range(ts) is used and the different timestamps uses
  # tifferent 'tzone' attributes, e.g. if some files where scanning during
  # daylight savings time and some not, we will get a warning saying:
  # "'tzone' attributes are inconsistent".  By doing the below, we avoid
  # this warning (which confuses users).
  ts <- sort(ts);
  ts <- ts[c(1,n)];
  ok <- is.finite(ts);
  if (any(ok)) {
    # range() gives strange values?!? Instead:
    ts[ok] <- format(ts[ok], "%Y-%m-%d %H:%M:%S");
  }
  s <- c(s, sprintf("Time period: %s -- %s", ts[1], ts[2]));
  s <- c(s, sprintf("Total file size: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("getTimestamps", "CnagCfhSet", function(this, ..., force=FALSE) {
  ts <- this$.timestamps;

  if (force || is.null(ts)) {
    # Get CFH header dates
    ts <- lapply(this, FUN=getTimestamp);
    ts <- do.call(c, args=ts);
    this$.timestamps <- ts;
  }

  ts;
})


setMethodS3("getIdentifier", "CnagCfhSet", function(this, ..., force=FALSE) {
  identifier <- this$.identifier;
  if (force || is.null(identifier)) {
    identifier <- NextMethod("getIdentifier");
    if (is.null(identifier)) {
      identifiers <- lapply(this, FUN=getIdentifier);
      identifier <- getChecksum(identifiers);
    }
    this$.identifier <- identifier;
  }
  identifier;
}, private=TRUE)




###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this CFH set"
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
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "CnagCfhSet", function(this, ...) {
  getCdf(this$files[[1]], ...);
})


###########################################################################/**
# @RdocMethod setCdf
#
# @title "Sets the CDF structure for this CFH set"
#
# \description{
#  @get "title".  This structures is assigned to all CFH files in the set.
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \seealso{
#   @seemethod "getCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("setCdf", "CnagCfhSet", function(this, cdf, verbose=FALSE, ...) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Setting CDF for CFH set");
  verbose && print(verbose, cdf);

  # Nothing to do?
#  oldCdf <- getCdf(this);
#  if (equals(cdf, oldCdf))
#    return(invisible(this));

  # Set the CDF for all CFH files
  verbose && enter(verbose, "Setting CDF for each CFH file");
  lapply(this, FUN=setCdf, cdf, ...);
  verbose && exit(verbose);

  # Have to clear the cache
  verbose && enter(verbose, "Clearing data-set cache");
  clearCache(this);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(this);
})


setMethodS3("findByName", "CnagCfhSet", function(static, ..., paths="cnagData(|,.*)/") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'paths':
  if (is.null(paths)) {
    paths <- eval(formals(findByName.CnagCfhSet)[["paths"]]);
  }

  NextMethod("findByName", paths=paths);
}, static=TRUE, protected=TRUE)


setMethodS3("byName", "CnagCfhSet", function(static, name, tags=NULL, chipType, paths=NULL, ...) {
  suppressWarnings({
    path <- findByName(static, name, tags=tags, chipType=chipType, paths=paths, ...);
  })
  if (is.null(path)) {
    path <- file.path(paste(c(name, tags), collapse=","), chipType);
    throw("Cannot create ", class(static)[1], ".  No such directory: ", path);
  }

  suppressWarnings({
    byPath(static, path=path, ...);
  })
}, static=TRUE)


setMethodS3("byPath", "CnagCfhSet", function(static, path="rawData/", pattern="[.](c|C)(f|F)(h|H)$", checkChipType=TRUE, ..., onDuplicates=c("keep", "exclude", "error"), fileClass="CnagCfhFile", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'onDuplicates':
  onDuplicates <- match.arg(onDuplicates);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Defining ", class(static)[1], " from files");

  ## Don't explicitly pass the first argument after 'static', otherwise
  ## it (here argument 'path') may be part of '...' as well. /HB 2013-07-28
  this <- NextMethod("byPath", pattern=pattern, fileClass=fileClass, verbose=less(verbose));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Handle duplicates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (onDuplicates %in% c("exclude", "error")) {
    dups <- isDuplicated(this);
    ndups <- sum(dups);
    if (ndups > 0) {
      dupsStr <- paste(getNames(this)[dups], collapse=", ");
      if (onDuplicates == "error") {
        msg <- paste("Detected ", ndups, " duplicated CFH files (same datestamp): ", dupsStr, sep="");
        verbose && cat(verbose, "ERROR: ", msg);
        throw(msg);
      } else if (onDuplicates == "exclude") {
        this <- extract(this, !dups, onDuplicates="error");
        msg <- paste("Excluding ", ndups, " duplicated CFH files (same datestamp): ", dupsStr, sep="");
        verbose && cat(verbose, "WARNING: ", msg);
        warning(msg);
      }
    }
  }

  if (length(this) > 0) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Scan all CFH files for possible chip types
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Chip type according to the directory structure
    path <- getPath(this);
    chipType <- basename(path);
    verbose && cat(verbose, "The chip type according to the directory is: ",
                                                                chipType);

    # Let the directory name specify the chip type?
    if (checkChipType) {
      verbose && enter(verbose, "Scanning CFH set for used chip types");
      # This takes time if the CFH files are ASCII files.
      chipTypes <- sapply(this, FUN=function(file) {
        readCfhHeader(getPathname(file))$chipType;
      })
      tChipTypes <- table(chipTypes);
      verbose && print(verbose, tChipTypes);
      nbrOfChipTypes <- length(tChipTypes);
      verbose && exit(verbose);

      if (nbrOfChipTypes > 1) {
        verbose && cat(verbose, "Detected ", nbrOfChipTypes,
                                        " different chip types in CFH set: ",
                                    paste(names(tChipTypes), collapse=", "));
      } else {
        verbose && cat(verbose, "All CFH files use the same chip type: ",
                                                          names(tChipTypes));
      }

      # If chip type is taken from CFH headers and there are more than
      # one chip type in the set, then it is an error.
      if (nbrOfChipTypes > 1) {
        throw("Detected ", nbrOfChipTypes, " different chip types in CFH set. Use argument 'checkChipType=FALSE' to let the directory name of the CFH set specify the chip type instead: ", paste(names(tChipTypes), collapse=", "));
      }

      # Validate that the directory name matches the chip type
      if (!identical(names(tChipTypes), chipType)) {
        throw("Invalid name of directory containing CFH files. The name of the directory (", chipType, ") must be the same as the chip type used for the CFH files (", names(tChipTypes), "): ", path);
      }
    } else {
      verbose && cat(verbose, "Since 'checkChipType=FALSE', then the chip type specified by the directory name is used: ", chipType);
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Use the same CDF object for all CFH files.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Updating the CDF for all files");
    verbose && cat(verbose, "Chip type: ", chipType);
    cdf <- AffymetrixCdfFile$byChipType(chipType);
    setCdf(this, cdf, .checkArgs=FALSE);
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Scan for SAF files and apply them
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Scanning for and applying sample annotation files");
    sas <- SampleAnnotationSet$loadAll(verbose=less(verbose));
    if (length(sas) == 0L) {
      verbose && cat(verbose, "No sample annotation files found.");
    } else {
      verbose && print(verbose, sas);
      setAttributesBy(this, sas);
    }
    # Store the SAFs for now.
    this$.sas <- sas;
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  this;
}, static=TRUE, protected=TRUE)




###########################################################################/**
# @RdocMethod as.CnagCfhSet
# @alias as.CnagCfhSet.list
# @alias as.CnagCfhSet.default
#
# @title "Coerce an object to an CnagCfhSet object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Other arguments passed to @see "base::list.files".}
# }
#
# \value{
#   Returns an @see "CnagCfhSet" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.CnagCfhSet", "CnagCfhSet", function(object, ...) {
  object;
})

setMethodS3("as.CnagCfhSet", "list", function(object, ...) {
  CnagCfhSet(object, ...);
})

setMethodS3("as.CnagCfhSet", "default", function(object, ...) {
  throw("Cannot coerce object to an CnagCfhSet object: ", mode(object));
})



###########################################################################/**
# @RdocMethod isDuplicated
#
# @title "Identifies duplicated CFH files"
#
# \description{
#   @get "title" by comparing the timestamps in the CFH headers.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @logical @vector of length equal to the number of files
#   in the set.
#   An element with value @TRUE indicates that the corresponding CFH file
#   has the same time stamp as another preceeding CFH file.
# }
#
# \seealso{
#   Internally @see "base::duplicated" is used to compare timestamps.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isDuplicated", "CnagCfhSet", function(this, ...) {
  # Get the CFH header timestamp for all files
  timestamps <- getTimestamps(this);

  dups <- duplicated(timestamps);
  names(dups) <- getNames(this);

  dups;
}, protected=TRUE)



setMethodS3("getData", "CnagCfhSet", function(this, indices=NULL, fields=c("x", "y", "intensities", "stdvs", "pixels"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indices':
  nbrOfCells <- nbrOfCells(getCdf(this));
  if (!is.null(indices)) {
    indices <- Arguments$getIndices(indices, max=nbrOfCells);
    nbrOfCells <- length(indices);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  nbrOfArrays <- length(this);
  verbose && enter(verbose, "Getting cell data for ", nbrOfArrays, " arrays.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating the return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating the return structure");
  nbrOfFields <- length(fields);
  res <- vector("list", nbrOfFields);
  names(res) <- fields;
  for (field in fields) {
    if (field %in% c("x", "y", "pixels")) {
      naValue <- as.integer(NA);
    } else {
      naValue <- as.double(NA);
    }
    res[[field]] <- matrix(naValue, nrow=nbrOfCells, ncol=nbrOfArrays);
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading cell signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving data from ", nbrOfArrays, " arrays");
  for (kk in seq_len(nbrOfArrays)) {
    verbose && enter(verbose, "Array #", kk, " of ", nbrOfArrays);
    dataFile <- this$files[[kk]];
    value <- getData(dataFile, indices=indices, fields=fields, verbose=less(verbose));
    for (field in fields) {
      res[[field]][,kk] <- value[[field]];
      value[[field]] <- NULL;
    }
    # Not needed anymore
    value <- NULL; gc();
    verbose && exit(verbose);
  }
  verbose && exit(verbose);


  verbose && exit(verbose);

  res;
}) # getData()



setMethodS3("readUnits", "CnagCfhSet", function(this, units=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "readCelUnits() of CnagCfhSet");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Generating hashcode key for cache");
  key <- list(method="readUnits", class=class(this)[1]);
  if (is.list(units)) {
    key <- c(key, units=names(units), ...);
  } else {
    key <- c(key, units=units, ...);
  }
  id <- getChecksum(key);
  verbose && exit(verbose);
  if (!force) {
    verbose && enter(verbose, "Trying to obtain cached data");
    res <- this$.readUnitsCache[[id]];
    verbose && exit(verbose);
    if (!is.null(res)) {
      verbose && cat(verbose, "readUnits(): Returning cached data");
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the pathnames of all CFH files
  pathnames <- getPathnames(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data from file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling readCelUnits() for ",
                                              length(pathnames), " files");
  if (is.list(units)) {
    res <- .readCelUnits(pathnames, cdf=units, dropArrayDim=FALSE, ...);
  } else {
    # Always ask for CDF information from the CDF object!
    verbose && enter(verbose, "Retrieving CDF unit information");
    suppressWarnings({
      cdf <- readUnits(getCdf(this), units=units, ..., verbose=less(verbose));
    });
    verbose && str(verbose, cdf[1]);
    verbose && exit(verbose);
    verbose && enter(verbose, "Retrieving CFH units across samples");
    res <- .readCelUnits(pathnames, cdf=cdf, dropArrayDim=FALSE, ...);
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
    verbose && cat(verbose, "readUnits(): Updated cache");
  }

  verbose && exit(verbose);

  res;
})


###########################################################################/**
# @RdocMethod getAverageFile
#
# @title "Calculates the mean and the standard deviation of the cell signal (intensity, standard deviation etc.) across the CFH set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{The label of the calculated parameters.
#    If @NULL, a default name format \code{<prefix>-<mean>-<sd>} is used.}
#  \item{indices}{An @integer @vector specifying which cells to consider.
#    If \code{"remaining"}, only parameters for cells that have not been
#    are calculated.
#    If @NULL, all cells are used.}
#  \item{mean}{A @character of a @function specifying the function used
#    to calculate the average.}
#  \item{sd}{A @character of a @function specifying the function used
#    to calculate the standard deviation.}
#  \item{na.rm}{If @TRUE, @NAs are excluded before, otherwise not.}
#  \item{...}{Not used.}
#  \item{cellsPerChunk}{A @integer specifying the total number of cells
#    (across arrays) read into memory per chunk.}
#  \item{moreCells}{A @double scalar indicating if more or less cells
#    should be used per chunk.}
#  \item{force}{If @TRUE, parameters for cells already calculated are
#    recalculated, otherwise not.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns an @see "CnagCfhSet" of the same class as the CFH set
#   averaged.
# }
#
# \details{
#   The parameter estimates are stored as a CFH file of the same class as
#   the data files in the set.  The CFH file is named \code{<name>.cel}
#   and placed in the directory of the set.
#   Currently there is no specific data class for this file, but the average
#   cell signals are stored as "intensities", the standard deviation of the
#   cell signals as "stddevs", and the number of data points used for each
#   estimate is stored as "pixels".
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAverageFile", "CnagCfhSet", function(this, name=NULL, prefix="average", indices="remaining", field=c("intensities", "stdvs"), mean=c("median", "mean"), sd=c("mad", "sd"), na.rm=TRUE, g=NULL, h=NULL, ..., cellsPerChunk=moreCells*10^7/length(this), moreCells=1, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'field':
  field <- match.arg(field);

  # Argument 'mean':
  if (is.character(mean)) {
    mean <- match.arg(mean);
    meanName <- mean;
    if (mean == "mean") {
      mean <- base::rowMeans;
    } else if (mean == "median") {
      mean <- rowMedians;
    }
  } else if (is.function(mean)) {
    meanName <- "customMean";
  } else {
    throw("Argument 'mean' must be either a character or a function: ", mode(mean));
  }

  # Argument 'sd':
  if (is.character(sd)) {
    sd <- match.arg(sd);
    sdName <- sd;
    if (sd == "sd") {
      sd <- rowSds;
    } else if (sd == "mad") {
      sd <- rowMads;
    }
  } else if (is.function(sd)) {
    sdName <- "customSd";
  } else {
    throw("Argument 'sd' must be either a character or a function: ",
                                                           mode(sd));
  }

  # Argument 'name':
  if (is.null(name)) {
    key <- list(method="getAverageFile", class=class(this)[1],
                arrays=sort(getNames(this)), mean=meanName, sd=sdName);
    # assign mean and sd to an empty environment so that digest() doesn't
    # pick up any "promised" objects from the original environment.
    # A bit ad hoc, but it works for now. /2007-01-03
    key <- lapply(key, FUN=function(x) {
      if (is.function(x))
        environment(x) <- emptyenv();
      x;
    })
    id <- getChecksum(key);
    name <- sprintf("%s-%s-%s-%s,%s", prefix, field, meanName, sdName, id);
  }

  # Argument 'indices':
  df <- as.list(this)[[1]];
  nbrOfCells <- getHeader(df)$total;
  if (force) {
    if (identical(indices, "remaining")) {
      indices <- NULL;
    }
  }

  if (is.null(indices)) {
    indices <- 1:nbrOfCells;
  } else if (identical(indices, "remaining")) {
  } else {
    indices <- Arguments$getIndices(indices, max=nbrOfCells);
  }

  # Argument 'cellsPerChunk':
  cellsPerChunk <- Arguments$getInteger(cellsPerChunk, range=c(1,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving average cell signals across ", length(this), " arrays");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create CFH file to store the average array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a private filename (with a dot prefix) to make sure it is not
  # identified as a regular CFH file when the directory is scanned for files.
  filename <- sprintf(".%s.CFH", name);
  if (is.null(this$.averageFiles))
    this$.averageFiles <- list();
  res <- this$.averageFiles[[filename]];

  if (is.null(res)) {
    verbose && enter(verbose, "Obtaining an (existing or new) result file");

    # Searching for the output file in multiple directories
    path <- getPath(this);
    paths <- c(path);

    # Drop tags from root path?
    if (getOption(aromaSettings, "devel/dropRootPathTags", TRUE)) {
      path <- dropRootPathTags(path, depth=2, verbose=less(verbose, 5));
      paths <- c(paths, path);
      paths <- unique(paths);
    }

    verbose && cat(verbose, "Paths:");
    verbose && print(verbose, paths);
    verbose && cat(verbose, "Filename: ", filename);

    pathname <- NULL;
    for (kk in seq_along(paths)) {
      path <- paths[kk];
      verbose && enter(verbose, sprintf("Searching path #%d of %d", kk, length(paths)));

      verbose && cat(verbose, "Path: ", path);
      pathnameT <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
      verbose && cat(verbose, "Pathname: ", pathnameT);
      if (isFile(pathnameT)) {
        pathname <- pathnameT;
        verbose && cat(verbose, "Found an existing file.");
        verbose && exit(verbose);
        break;
      }

      verbose && exit(verbose);
    } # for (kk ...)
    verbose && cat(verbose, "Located pathname: ", pathname);

    verbose && exit(verbose);

    if (isFile(pathname)) {
      verbose && enter(verbose, "Loading existing data file");
      verbose && cat(verbose, "Pathname: ", pathname);
      res <- newInstance(df, pathname);
      verbose && exit(verbose);
    } else {
      verbose && enter(verbose, "Creating CFH file to store average signals");
      path <- paths[length(paths)];

      verbose && cat(verbose, "Path: ", path);
      verbose && cat(verbose, "Filename: ", filename);

      res <- createFrom(df, filename=filename, path=path,
                        methods="create", clear=TRUE, verbose=less(verbose));

      verbose && exit(verbose);
    } # if (isFile(pathname))

    this$.averageFiles[[filename]] <- res;
  } # if (is.null(res))

  verbose && print(verbose, res);

  pathname <- getPathname(res);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify which indices to use
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(indices, "remaining")) {
    pixels <- .readCel(pathname, readIntensities=FALSE, readStdvs=FALSE,
                      readPixels=TRUE)$pixels;
    indices <- which(pixels == 0);
    # Not needed anymore
    pixels <- NULL; # Not needed anymore.
  }

  nbrOfIndices <- length(indices);

  # Nothing more to do?
  if (nbrOfIndices == 0)
    return(res);

  verbose && cat(verbose, "Number of cells to be updated: ", nbrOfIndices);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate the mean and standard deviation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Since we might want to do this robustly, but also because we want to
  # estimate the standard deviation, for each cell we need all data across
  # arrays at once.  In order to this efficiently, we do this in chunks
  idxs <- 1:nbrOfIndices;
  head <- 1:cellsPerChunk;
  nbrOfChunks <- ceiling(nbrOfIndices / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  # Get the pathnames of all CFH files to average
  pathnames <- lapply(this, FUN=getPathname);
  pathnames <- unlist(pathnames, use.names=FALSE);
  nbrOfArrays <- length(pathnames);

  if (!na.rm)
    n <- rep(nbrOfArrays, length=cellsPerChunk);
  count <- 1;
  while (length(idxs) > 0) {
    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < cellsPerChunk) {
      head <- 1:length(idxs);
      if (!na.rm)
        n <- rep(nbrOfArrays, length=length(idxs));
    }

    # The indices to be used in this chunk
    ii <- idxs[head];
    verbose && cat(verbose, "Chunk size: ", length(ii));

    verbose && enter(verbose, "Reading data");
#    X <- readCelIntensities(pathnames, indices=indices[ii]);
    readIntensities <- (field == "intensities");
    readStdvs <- (field == "stdvs");
    # TODO: Ideally, affxparser::readCel() should support
    # multiple filenames turning every data fields into a
    # matrix. /HB 2007-01-07
    naValue <- as.double(NA);
    X <- matrix(naValue, nrow=length(ii), ncol=nbrOfArrays);
    for (kk in seq_len(nbrOfArrays)) {
      X[,kk] <- .readCel(filename = pathnames[kk],
                        indices = indices[ii],
                        readIntensities = readIntensities,
                        readHeader = FALSE,
                        readStdvs = readStdvs,
                        readPixels = FALSE,
                        readXY = FALSE,
                        readOutliers = FALSE,
                        readMasked = FALSE,
                        ...,
                        verbose = (verbose - 1))[[field]];
    }
    verbose && exit(verbose);

    if (!is.null(g)) {
      verbose && enter(verbose, "Transforming data using y = g(x)");
      X <- g(X);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Estimating averages and standard deviations");
    if (na.rm)
      n <- base::apply(X, MARGIN=1, FUN=function(x) { sum(!is.na(x)) });

    # Calculate the mean signal
    mu <- mean(X, na.rm=na.rm);          # Special mean()!
    # Calculate the standard deviation of the signals
    sigma <- sd(X, mean=mu, na.rm=na.rm);   # Special sd()!

    verbose && exit(verbose);

    if (!is.null(h)) {
      verbose && enter(verbose, "Back-transforming estimates using x = h(y)");
      mu <- h(mu);
      sigma <- h(sigma);
      verbose && exit(verbose);
    }

    # Write estimates to result file
    verbose && enter(verbose, "Writing estimates");
    .updateCel(pathname, indices=indices[ii], intensities=mu, stdvs=sigma, pixels=n);
    verbose && exit(verbose);

    # Not needed anymore
    mu <- sigma <- NULL;

    # Next chunk...
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    gc();
    verbose && exit(verbose);
  } # while()

  verbose && exit(verbose);

  res;
})


setMethodS3("getAverage", "CnagCfhSet", function(this, ...) {
  getAverageFile(this, ...);
})

setMethodS3("getAverageLog", "CnagCfhSet", function(this, ...) {
  getAverageFile(this, g=log2, h=function(x) 2^x, ...);
})

setMethodS3("getAverageAsinh", "CnagCfhSet", function(this, ...) {
  getAverageFile(this, g=asinh, h=sinh, ...);
})



setMethodS3("range", "CnagCfhSet", function(this, ...) {
  range(unlist(lapply(this, FUN=range, ...), use.names=FALSE));
}, protected=TRUE)


setMethodS3("getDefaultFullName", "CnagCfhSet", function(this, parent=1L, ...) {
  NextMethod("getDefaultFullName", parent=parent);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-11-20
# o CLEANUP: Deprecated "[" and "[[", because they should be used to
#   subset files and not units.
# 2011-03-03
# o GENERALIZATION: Now CnagCfhSet locates sample annotation files and
#   sets the attributes following the new aroma search convention.
# 2011-02-28
# o Now getAverageFile() first tries to locate an existing result file
#   in multiple root paths.  If not found, it creates a new one.
# 2011-02-24
# o GENERALIZATION: Now getAverageFile() for CnagCfhSet drops tags
#   from the output root path (if 'devel/dropRootPathTags' setting is TRUE).
# o Expanded the searched root paths to be cnagData(|,.*)/.
# 2009-08-12
# o Now findByName() of CnagCfhSet utilizes ditto of AffymetrixCelSet,
#   because its code was identical to the latter.
# 2008-07-21
# o Now findByName() assert that the data set name is not empty.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: getData() and getAverageFile().
# 2008-05-08
# o If paths=NULL in findByName(), it becomes the default argument value.
# 2007-08-10
# o Now getAverageFile() utilizes Biobase::rowMedians(), if available.
# 2007-04-06
# o Created from AffymetrixCelSet.R
############################################################################
