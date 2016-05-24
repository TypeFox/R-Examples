###########################################################################/**
# @RdocClass AffymetrixCelSet
#
# @title "The AffymetrixCelSet class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixCelSet object represents a set of Affymetrix CEL files
#  with \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "AffymetrixCelFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \examples{\dontrun{
#   @include "../incl/AffymetrixCelSet.Rex"
# }}
#
# \seealso{
#   @see "AffymetrixCelFile".
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AffymetrixCelSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    reqFileClass <- "AffymetrixCelFile";
    lapply(files, FUN=function(df) {
      df <- Arguments$getInstanceOf(df, reqFileClass, .name="files");
    })
  } else if (inherits(files, "AffymetrixCelSet")) {
    return(as.AffymetrixCelSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  this <- extend(AffymetrixFileSet(files=files, ...), "AffymetrixCelSet",
    "cached:.intensities" = NULL,
    "cached:.intensitiesIdxs" = NULL,
    "cached:.readUnitsCache" = NULL,
    "cached:.getUnitIntensitiesCache" = NULL,
    "cached:.averageFiles" = list(),
    "cached:.timestamps" = NULL,
    "cached:.fileSize" = NULL
  );

  if (length(this$files) > 0) {
    # Make sure the set name is non-empty
    name <- getName(this);
    if (!nzchar(name)) {
      throw("An ", class(this)[1], " must have a name of at least length one: ", getPathname(this))
    }
  }

  this;
})


setMethodS3("clearCache", "AffymetrixCelSet", function(this, ...) {
  # Then for this object
  NextMethod("clearCache");

  this$.averageFiles <- list();

  if (length(this) > 0) {
    # Clear the cache for the CDF.
    cdf <- getCdf(this);
    clearCache(cdf);
  }

  invisible(this);
}, private=TRUE)


setMethodS3("clone", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Cloning Affymetrix CEL set");

  # Clone itself and the files.  The call below will clear the cache!
  object <- NextMethod("clone", clear=TRUE, verbose=less(verbose));
  clearCache(object);

  if (length(object) > 0) {
    # Clone the CDF (this will update the CDF of all file object)
    verbose && enter(verbose, "Cloning CDF");
    cdf <- getCdf(object);
    cdf <- clone(cdf);
    verbose && exit(verbose);
    verbose && enter(verbose, "Adding CDF to CEL set");
    setCdf(object, cdf, .checkArgs=FALSE);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  object;
}, protected=TRUE)



setMethodS3("append", "AffymetrixCelSet", function(this, other, ..., verbose=FALSE) {
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, class(this)[1]);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Appending CEL set");
  verbose && print(verbose, other);

  # Validate chip type
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  for (file in getFiles(other)) {
    oCdf <- getCdf(file);
    oChipType <- getChipType(oCdf);
    if (!identical(oChipType, chipType)) {
      throw("Argument 'other' contains a CEL file of different chip type: ",
                                                oChipType, " != ", chipType);
    }
  }


  # Append using append() method in super class
  # NOTE: Do not pass arguments, they haven't been modified
  # and do not need to be passed. Indeed, it will give an error
  # for unknown reasons, cf. aroma.affymetrix thread 'Append function'
  # on March 10, 2011. /HB 2011-03-10
  this <- NextMethod("append");

  # Set the same CDF for all CEL files
  verbose && enter(verbose, "Updating the CDF for all files");
  setCdf(this, cdf);
  verbose && exit(verbose);

  verbose && exit(verbose);

  this;
}, protected=TRUE)




###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the Affymetrix CEL set"
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
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "AffymetrixCelSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  s <- c(s, sprintf("Tags: %s", tags));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  n <- length(this);
  s <- c(s, sprintf("Number of arrays: %d", n));
  names <- getNames(this);
  s <- c(s, sprintf("Names: %s [%d]", hpaste(names), n));

  # Get CEL header timestamps?
  maxCount <- getOption(aromaSettings, "output/timestampsThreshold");
  if (n == 0) {
    s <- c(s, "Time period: NA");
  } else if (maxCount >= n) {
    ts <- getTimestamps(this);
    # Note: If ts <- range(ts) is used and the different timestamps uses
    # tifferent 'tzone' attributes, e.g. if some files where scanning during
    # daylight savings time and some not, we will get a warning saying:
    # "'tzone' attributes are inconsistent".  By doing the below, we avoid
    # this warning (which confuses users).
    ts <- sort(ts);
    ts <- ts[c(1,n)];
    ts <- format(ts, "%Y-%m-%d %H:%M:%S");  # range() gives strange values?!?
    s <- c(s, sprintf("Time period: %s -- %s", ts[1], ts[2]));
  } else {
    s <- c(s, sprintf("Time period: [not reported if more than %.0f arrays]", as.double(maxCount)));
  }
  s <- c(s, sprintf("Total file size: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("getTimestamps", "AffymetrixCelSet", function(this, ..., force=FALSE) {
  ts <- this$.timestamps;

  if (force || is.null(ts)) {
    # Get CEL header dates
    ts <- lapply(this, FUN=getTimestamp);
    ts <- do.call(c, args=ts);
    this$.timestamps <- ts;
  }

  ts;
})


setMethodS3("getIdentifier", "AffymetrixCelSet", function(this, ..., force=FALSE) {
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
# @RdocMethod getChipType
#
# @title "Gets the chip type for this CEL set"
#
# \description{
#  @get "title".
#  The chip type is inferred from the current CDF.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @character string.
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
setMethodS3("getChipType", "AffymetrixCelSet", function(this, ...) {
  if (length(this) == 0) {
    return(as.character(NA));
  }
  unf <- getUnitNamesFile(this);
  getChipType(unf, ...);
})

setMethodS3("getPlatform", "AffymetrixCelSet", function(this, ...) {
  "Affymetrix";
})


setMethodS3("getUnitNamesFile", "AffymetrixCelSet", function(this, ...) {
  getUnitNamesFile(getOneFile(this), ...);
})

setMethodS3("getUnitTypesFile", "AffymetrixCelSet", function(this, ...) {
  getUnitTypesFile(getOneFile(this), ...);
})



###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this CEL set"
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
setMethodS3("getCdf", "AffymetrixCelSet", function(this, ...) {
  getCdf(getOneFile(this), ...);
})


###########################################################################/**
# @RdocMethod setCdf
#
# @title "Sets the CDF structure for this CEL set"
#
# \description{
#  @get "title".  This structures is assigned to all CEL files in the set.
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#     May also be a @see "R.utils::Verbose" object.}
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
setMethodS3("setCdf", "AffymetrixCelSet", function(this, cdf, verbose=FALSE, ..., .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'cdf':
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

    # Assure that the CDF is compatible with the CEL file
    if (length(this) > 0L) {
      cf <- getOneFile(this);
      if (nbrOfCells(cdf) != nbrOfCells(cf)) {
        throw("Cannot set CDF. The specified CDF structure ('", getChipType(cdf), "') is not compatible with the chip type ('", getChipType(cf), "') of the CEL file. The number of cells do not match: ", nbrOfCells(cdf), " != ", nbrOfCells(cf));
      }
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Setting CDF for CEL set");
  verbose && print(verbose, cdf);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # ASCII CDFs are only allowed if explicited accepted in the rules.
  # To change, do:
  #  setOption(aroma.affymetrix, "rules$allowAsciiCdfs", TRUE)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allowAsciiCdfs <- getOption(aromaSettings, "rules/allowAsciiCdfs", FALSE);
  if (!allowAsciiCdfs) {
    # ASCII CDF are *not* allowed
    ff <- getFileFormat(cdf);
    if (regexpr("ASCII", ff) != -1) {
      throw("Cannot set CDF for data set. The given CDF is in ASCII format, which is protected against by default. It is much faster and more memory efficient and *highly* recommended to work with binary CDF files. Use affxparser::convertCdf() to convert a CDF into another format.  For advanced users (who truly know what they are doing) it is indeed possible to use ASCII CDFs (see website for details). Details on the CDF file: ", getPathname(cdf), " [", ff, "].");
    }
  }

  # Set the CDF for all CEL files
  verbose && enter(verbose, "Setting CDF for each CEL file");
  lapply(this, FUN=setCdf, cdf, .checkArgs=FALSE, ...);
  verbose && exit(verbose);

  # Have to clear the cache
  verbose && enter(verbose, "Clearing data-set cache");
  clearCache(this);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(this);
})


setMethodS3("findByName", "AffymetrixCelSet", function(static, ..., chipType=NULL, paths=c("rawData(|,.*)/", "probeData(|,.*)/")) {
  # Arguments 'chipType':`
  if (!is.null(chipType)) {
    chipType <- Arguments$getCharacter(chipType);
  }

  # Arguments 'paths':
  if (is.null(paths)) {
    paths <- eval(formals(findByName.AffymetrixCelSet)[["paths"]]);
  }

  # Call the "next" method
  NextMethod("findByName", subdirs=chipType, paths=paths);
}, static=TRUE, protected=TRUE)



setMethodS3("byName", "AffymetrixCelSet", function(static, name, tags=NULL, chipType=NULL, cdf=NULL, paths=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  if (!is.null(chipType)) {
    chipType <- Arguments$getCharacter(chipType);
  }

  # Argument 'cdf':
  if (!is.null(cdf)) {
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Setting up ", class(static)[1], " by name");

  verbose && cat(verbose, "Name: ", name);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=","));

  if (is.null(cdf) && is.null(chipType)) {
    throw("Either argument 'chipType' or argument 'cdf' must be specified.");
  }

#  if (!is.null(cdf) && !is.null(chipType)) {
#    throw("Only one of argument 'chipType' and argument 'cdf' can be specified.");
#  }

  # If argument 'cdf' is not specified, infer it from 'chipType'.
  if (is.null(cdf)) {
    cdf <- AffymetrixCdfFile$byChipType(chipType);
  }

  # The chiptype without tags
  chipTypeShort <- getChipType(cdf, fullname=FALSE);

  verbose && cat(verbose, "Chip type: ", chipTypeShort);

  suppressWarnings({
    paths <- findByName(static, name, tags=tags, chipType=chipTypeShort,
                       paths=paths, firstOnly=FALSE, ...);
  })
  if (is.null(paths)) {
    path <- file.path(paste(c(name, tags), collapse=","), chipTypeShort);
    throw("Cannot create ", class(static)[1], ".  No such directory: ", path);
  }

  verbose && cat(verbose, "Paths to possible data sets:");
  verbose && print(verbose, paths);

  # Record all exception
  exList <- list();

  res <- NULL;
  for (kk in seq_along(paths)) {
    path <- paths[kk];
    verbose && enter(verbose, sprintf("Trying path #%d of %d", kk, length(paths)));
    verbose && cat(verbose, "Path: ", path);

    tryCatch({
      suppressWarnings({
        res <- byPath(static, path=path, cdf=cdf, ..., verbose=verbose);
      });
    }, error = function(ex) {
      exList <<- append(exList, list(list(exception=ex, path=path)));

      verbose && cat(verbose, "Data set could not be setup for this path, because:");
      verbose && cat(verbose, ex$message);
    });

    if (!is.null(res)) {
      if (length(res) > 0) {
        verbose && cat(verbose, "Successful setup of data set.");
        verbose && exit(verbose);
        break;
      }
    }

    verbose && exit(verbose);
  } # for (kk ...)

  if (is.null(res)) {
    exMsgs <- sapply(exList, FUN=function(ex) {
      sprintf("%s (while trying '%s').",
                   ex$exception$message, ex$path);
    });
    exMsgs <- sprintf("(%d) %s", seq_along(exMsgs), exMsgs);
    exMsgs <- paste(exMsgs, collapse="  ");
    msg <- sprintf("Failed to setup a data set for any of %d data directories located. The following reasons were reported: %s", length(paths), exMsgs);
    verbose && cat(verbose, msg);
    throw(msg);
  }

  verbose && exit(verbose);

  res;
}, static=TRUE) # byName()



setMethodS3("update2", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating ", class(this)[1]);
  updateSampleAnnotationSet(this, ..., verbose=less(verbose, 1));
  verbose && exit(verbose);
}, protected=TRUE)


setMethodS3("updateSampleAnnotationSet", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scan for SAF files and apply them
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Scanning for and applying sample annotation files");

  # Nothing to do?
  if (length(this) > 0) {
    sas <- SampleAnnotationSet$loadAll(verbose=less(verbose));
    if (length(sas) == 0) {
      verbose && cat(verbose, "No sample annotation files found.");
    } else {
      verbose && print(verbose, sas);
      setAttributesBy(this, sas);
    }
    # Store the SAFs for now.
    this$.sas <- sas;
  } else {
    verbose && cat(verbose, "Empty data set. Nothing to do.");
  }

  verbose && exit(verbose);

  invisible(this);
}, private=TRUE)  # updateSampleAnnotationSet()


setMethodS3("byPath", "AffymetrixCelSet", function(static, path, cdf=NULL, pattern="[.](c|C)(e|E)(l|L)$", ..., checkChipType=is.null(cdf), onDuplicates=c("keep", "exclude", "error"), fileClass="AffymetrixCelFile", force=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'onDuplicates':
  onDuplicates <- match.arg(onDuplicates);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Defining ", class(static)[1], " from files");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Look for cached results (useful for extremely large data set)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="byPath", class=class(static)[1], path=path, pattern=pattern, cdf=cdf, checkChipType=checkChipType, ..., fileClass=fileClass);
  dirs <- c("aroma.affymetrix", "dataSets", class(static)[1]);
  res <- loadCache(key=key, dirs=dirs);
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "Found cached results");
    verbose && exit(verbose);
    return(res);
  }


  # Call the "next" method

  # WORKAROUND: For unknown reasons it is not possible to specify
  # 'path=path' below, because it is already passed implicitly by
  # NextMethod() and if done, then argument 'recursive' to byPath() for
  # GenericDataFileSet will also be assign the value of 'path', e.g.
  # try byPath(AffymetrixCelSet(), "path/to/").  This seems to be related
  # to R-devel thread 'Do *not* pass '...' to NextMethod() - it'll do it
  # for you; missing documentation, a bug or just me?' on Oct 16-18, 2012.
  args <- list("byPath", pattern=pattern, fileClass=fileClass, verbose=less(verbose));
#  set <- do.call(NextMethod, args);
  set <- NextMethod(generic="byPath", object=static, pattern=pattern, fileClass=fileClass, verbose=less(verbose));
  verbose && cat(verbose, "Retrieved files: ", length(set));

  if (length(set) > 0) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Handle duplicates
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (onDuplicates %in% c("exclude", "error")) {
      dups <- isDuplicated(set, verbose=less(verbose));
      ndups <- sum(dups);
      if (ndups > 0) {
        dupsStr <- paste(getNames(set)[dups], collapse=", ");
        if (onDuplicates == "error") {
          msg <- paste("Detected ", ndups, " duplicated CEL files (same datestamp): ", dupsStr, sep="");
          verbose && cat(verbose, "ERROR: ", msg);
          throw(msg);
        } else if (onDuplicates == "exclude") {
          set <- extract(set, !dups, onDuplicates="error");
          msg <- paste("Excluding ", ndups, " duplicated CEL files (same datestamp): ", dupsStr, sep="");
          verbose && cat(verbose, "WARNING: ", msg);
          warning(msg);
        }
      }
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Scan all CEL files for possible chip types
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Chip type according to the directory structure
    path <- getPath(set);
    chipType <- basename(path);
    verbose && cat(verbose,
                   "The chip type according to the path is: ", chipType);

    # Let the directory name specify the chip type?
    if (checkChipType) {
      verbose && enter(verbose, "Scanning CEL set for used chip types");
      # This takes time if the CEL files are ASCII files.
      chipTypes <- sapply(set, FUN=function(file) {
        .readCelHeader(getPathname(file))$chiptype;
      })
      tChipTypes <- table(chipTypes);
      verbose && print(verbose, tChipTypes);
      nbrOfChipTypes <- length(tChipTypes);
      verbose && exit(verbose);

      if (nbrOfChipTypes > 1) {
        verbose && cat(verbose, "Detected ", nbrOfChipTypes,
                                        " different chip types in CEL set: ",
                                    paste(names(tChipTypes), collapse=", "));
      } else {
        verbose && cat(verbose, "All CEL files use the same chip type: ",
                                                          names(tChipTypes));
      }

      # If chip type is taken from CEL headers and there are more than
      # one chip type in the set, then it is an error.
      if (nbrOfChipTypes > 1) {
        throw("Detected ", nbrOfChipTypes, " different chip types in CEL set. Use argument 'checkChipType=FALSE' to let the directory name of the CEL set specify the chip type instead: ", paste(names(tChipTypes), collapse=", "));
      }

      # Validate that the directory name matches the chip type
      tChipTypesShort <- gsub(",.*", "", names(tChipTypes));
      if (!identical(tChipTypesShort, chipType)) {
        throw("Invalid name of directory containing CEL files. The name of the directory (", chipType, ") must be the same as the chip type used for the CEL files (", names(tChipTypes), ") unless using argument 'checkChipType=FALSE': ", path);
      }
    } else {
      verbose && cat(verbose, "Since 'checkChipType=FALSE', then the chip type specified by the directory name is used: ", chipType);
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Use the same CDF object for all CEL files.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for chip type '", chipType, "' inferred from path");
      cf <- getOneFile(set);
      nbrOfCells <- nbrOfCells(cf);
      cdf <- AffymetrixCdfFile$byChipType(chipType, nbrOfCells=nbrOfCells);
      verbose && exit(verbose);

      verbose && enter(verbose, "Check compatibility with 1st CEL file");
      verbose && cat(verbose, "Chip type: ", chipType);
      if (nbrOfCells(cdf) != nbrOfCells(cf)) {
        cdf <- getCdf(cf);
        chipType <- getChipType(cdf);
        verbose && cat(verbose, "Chip type (updated): ", chipType);
      }
      verbose && exit(verbose);
    } else {
      verbose && cat(verbose, "Using prespecified CDF: ",
                     getChipType(cdf, fullname=TRUE));
    }
  }

  verbose && enter(verbose, "Updating the CDF for all files");
  if (!is.null(cdf)) {
    setCdf(set, cdf);
  }
  verbose && exit(verbose);

  # Let the new CEL set update itself
  update2(set, verbose=less(verbose, 1));

  # Save to file cache
  saveCache(set, key=key, dirs=dirs);

  verbose && exit(verbose);

  set;
}, static=TRUE, protected=TRUE) # byPath()




###########################################################################/**
# @RdocMethod as.AffymetrixCelSet
# @alias as.AffymetrixCelSet.list
# @alias as.AffymetrixCelSet.default
#
# @title "Coerce an object to an AffymetrixCelSet object"
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
#   Returns an @see "AffymetrixCelSet" object.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixCelSet", "AffymetrixCelSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixCelSet", "list", function(object, ...) {
  AffymetrixCelSet(object, ...);
})

setMethodS3("as.AffymetrixCelSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixCelSet object: ", mode(object));
})



###########################################################################/**
# @RdocMethod isDuplicated
#
# @title "Identifies duplicated CEL files"
#
# \description{
#   @get "title" by comparing the timestamps in the CEL headers.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a @logical @vector of length equal to the number of files
#   in the set.
#   An element with value @TRUE indicates that the corresponding CEL file
#   has the same time stamp as another preceeding CEL file.
# }
#
# @author "HB"
#
# \seealso{
#   Internally @see "base::duplicated" is used to compare timestamps.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isDuplicated", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Scanning for duplicated data files");

  verbose && enter(verbose, "Reading the timestamp of all files");
  # Get the CEL header timestamp for all files
  timestamps <- getTimestamps(this);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying duplicated timestamps");
  dups <- duplicated(timestamps);
  names(dups) <- getNames(this);
  verbose && exit(verbose);

  if (verbose) {
    if (any(dups)) {
      cat(verbose, "Duplicated files:", paste(names(dups), collapse=", "));
    } else {
      cat(verbose, "Duplicated files: none.");
    }
  }

  verbose && exit(verbose);

  dups;
}, protected=TRUE)



setMethodS3("getData", "AffymetrixCelSet", function(this, indices=NULL, fields=c("x", "y", "intensities", "stdvs", "pixels"), ..., verbose=FALSE) {
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
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


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


###########################################################################/**
# @RdocMethod getIntensities
#
# @title "Gets cell intensities from a set of cells and a set of arrays"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Passed to @seemethod "getData".}
# }
#
# \value{
#   Returns a @numeric \eqn{NxK} matrix, where \eqn{N} is the number of
#   cells read, and \eqn{K} is the number of arrays in the data set.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getIntensities", "AffymetrixCelSet", function(this, ...) {
  getData(this, ..., fields="intensities")$intensities;
}) # getIntensities()



###########################################################################/**
# @RdocMethod getUnitIntensities
#
# @title "Gets cell signals for a subset of units and a subset of arrays"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{An @integer index @vector specifying units to be read.
#    If @NULL, all units are read.}
#  \item{...}{Arguments passed to the low-level function for reading
#    CEL units, i.e. @see "affxparser::readCelUnits".
#    If \code{units} is not already a CDF @list structure, these arguments
#    are also passed to \code{readUnits()} of @see "AffymetrixCdfFile".}
#  \item{force}{If @TRUE, cached values are ignored.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a named @list structure.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getUnitIntensities", "AffymetrixCelSet", function(this, units=NULL, ..., force=FALSE, cache=!is.null(units), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="getUnitIntensities", class=class(this)[1], units=units, ...);
  id <- getChecksum(key);
  res <- this$.getUnitIntensitiesCache[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "getUnitIntensitiesCache(): Returning cached data");
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the pathnames of all CEL files
  pathnames <- unlist(lapply(this, FUN=getPathname), use.names=FALSE);

  # Is 'units' a (pre-created) CDF structure?
  if (is.list(units)) {
    cdfUnits <- units;
  } else {
    # Always ask for CDF information from the CDF object!
    cdf <- getCdf(this);
    cdfUnits <- getCellIndices(cdf, units=units, ...);
  }

  res <- .readCelUnits(pathnames, cdf=cdfUnits, readStdvs=FALSE,
                              readPixels=FALSE, dropArrayDim=FALSE, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    this$.getUnitIntensitiesCache <- list();
    this$.getUnitIntensitiesCache[[id]] <- res;
    verbose && cat(verbose, "readUnits(): Updated cache");
  }

  res;
})



setMethodS3("readUnits", "AffymetrixCelSet", function(this, units=NULL, ..., force=FALSE, cache=!is.null(units), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "readCelUnits() of AffymetrixCelSet");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Generating hashcode key for cache");
  key <- list(method="readUnits", class=class(this)[1]);
  if (is.list(units)) {
    # BUG FIX: Cannot stratify only on unit names, e.g. then 'stratifyBy'
    # will not make a difference! 2008-03-11
#    key <- c(key, units=names(units), ...);
#   ...will also not work (stratifyBy="pm" and "mm" gives the same)!
#    key <- c(key, unitNames=names(units), cdfUnitSize=object.size(units), ...);
    # This will allow us to pass pre-digested object.
    unitsHashcode <- attr(units, "hashcode");
    if (is.null(unitsHashcode)) {
      unitsHashcode <- getChecksum(units);
      attr(units, "hashcode") <- unitsHashcode;
    }
    key <- c(key, unitsHashcode=unitsHashcode, ...);
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
  # Get the pathnames of all CEL files
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
    cdf <- getCdf(this);
    cdfList <- getCellIndices(cdf, units=units, ..., verbose=less(verbose));
#    suppressWarnings({
#      cdfList <- readUnits(cdf, units=units, ..., verbose=less(verbose));
#    });
    verbose && str(verbose, cdfList[1]);
    verbose && exit(verbose);
    verbose && enter(verbose, "Retrieving CEL units across samples");
    res <- .readCelUnits(pathnames, cdf=cdfList, dropArrayDim=FALSE, ...);
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


setMethodS3("getAverage", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, ...);
})

setMethodS3("getAverageLog", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, g=log2, h=function(x) 2^x, ...);
})

setMethodS3("getAverageAsinh", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, g=asinh, h=sinh, ...);
})



setMethodS3("range", "AffymetrixCelSet", function(this, ...) {
  range(unlist(lapply(this, FUN=range, ...), use.names=FALSE));
}, protected=TRUE)



setMethodS3("applyToUnitIntensities", "AffymetrixCelSet", function(this, units=NULL, FUN, stratifyBy="pm", verbose=FALSE, ...) {
  y <- getUnitIntensities(this, units=units, stratifyBy=stratifyBy, ...);

  y <- lapply(y, FUN=function(unit) {
    groups <- lapply(unit, FUN=function(group) {
      FUN(group[[1]], ...)
    })
    groups;
  })
  y;
}, private=TRUE)


setMethodS3("getUnitGroupCellMap", "AffymetrixCelSet", function(this, ...) {
  getUnitGroupCellMap(getOneFile(this), ...);
})


############################################################################
# HISTORY:
# 2013-03-22
# o Obscufated the error message on ASCII CDFs on how to change settings
#   to allow the processing of such, because it shouldn't be easy to
#   do that.
# 2012-11-28
# o CLARIFICATION: Now the cache path used by byPath() for AffymetrixCelSet
#   includes "dataSets" and <staticClassName>.
# 2012-11-20
# o CLEANUP: Deprecated "[" and "[[", because they should be used to
#   subset files and not units.
# 2012-09-12
# o MEMORY: Now getUnitIntensities() and readUnits() for AffymetrixCelSet
#   no longer cache the results in memory if argument 'units' is NULL.
# o BUG FIX/ROBUSTNESS: Although getUnitIntensities() for AffymetrixCelSet
#   would retrieve the CDF tree structure and pass it to affxparser's
#   readCelIntensities(), it was done in such a way that the latter would
#   still need to find and read the CDF.  This meant that it would not
#   necessarily retrieve the same CDF as the AffymetrixCelSet used.
# 2011-08-16
# o ROBUSTNESS: If static byName() for AffymetrixCelSet fails to setup
#   a data set, it now reports the error message for each data set
#   directory it tried.
# 2011-03-03
# o GENERALIZATION: Now AffymetrixCelSet locates sample annotation files
#   and sets the attributes following the new aroma search convention.
# 2011-03-01
# o Now byName() handles if byPath() throws an exception.
# o BUG FIX: The recent updates of byName() forgot to pass argument 'cdf'
#   to the byPath() call.
# 2011-02-25
# o Now byName() for AffymetrixCelSet tries to setup all data set
#   directories matching the query and not just the first one.
#   It also requires that the data sets found by byPath() are non-empty,
#   otherwise they are skipped.
# 2011-02-24
# o Expanded the searched root paths to be rawData(|,.*)/ and
#   probeData(|,.*)/.
# 2011-02-04
# o BUG FIX: as.character() for an empty AffymetrixCelSet would throw
#   exception "Argument 'x' is of length 1 although the range ([0,0])
#   implies that is should be empty.", while trying to infer the chip type.
# 2010-07-02
# o Now AffymetrixCelSet$byName(..., chipType="GenomeWideSNP_6,Full") will
#   work (before chiptypes with tags would give an error).  This is now
#   done by first locating the CDF for the chip type (with tags).
# 2009-12-03
# o BUG FIX: The test for allowing ASCII CDFs or not in setCdf() of an
#   AffymetrixCelSet was only applied if getOption(aromaSettings,
#   "rules/allowAsciiCdfs") was TRUE (should be FALSE).
# 2009-08-12
# o Now findByName() of AffymetrixCelSet calls ditto of
#   GenericDataFileSet of the R.filesets packages.
# 2009-07-08
# o Added getUnitTypesFile() for AffymetrixCelSet.
# 2009-05-23
# o Now the chip type validation of fromFiles() for AffymetrixCelSet
#   is aware of tags in the chip type of the CEL files. This may happen
#   if custom CDFs are used are their full chip types are stored in the
#   CEL files, e.g. Hs_PromPR_v02,Harvard,ROIs,unique.
# 2008-12-18
# o BUG FIX: getUnitIntensities() of AffymetrixCelSet would drop the array
#   dimension if only one array was read.
# 2008-07-21
# o Now findByName() assert that the data set name is not empty.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: getData(), getAverageFile().
# 2008-07-14
# o Added protected update().
# o Added private updateSampleAnnotationSet().
# 2008-06-22
# o BUG FIX: Added missing getPlatform().
# 2008-05-31
# o BUG FIX: readUnits() would throw 'Error in readCelUnits(pathnames, cdf
#   = cdf, ...) : No CDF file for chip type found: GenomeWideSNP_6', if
#   the CDF was set to GenomeWideSNP_6,Full.CDF and no GenomeWideSNP_6.cdf
#   file was found.  This was because readUnits() retrieve (x,y) information
#   units, which requires affxparser::readCelUnits() to locate the CDF to
#   infer the number of probe column in order to map (x,y) to cell indices.
#   Now readUnits() get the cell indices directly using getCellIndices().
#   Thanks Yue Hu for noticing this problem.
# 2008-05-17
# o Added getChipType().
# 2008-05-09
# o Now AffymetrixCelSet inherits from AromaMicroarrayDataSet.
# 2008-05-08
# o Now we are using foo(static, ...) instead of static$foo(...).
# o If paths=NULL in findByName(), it becomes the default argument value.
# o BUG FIX: fromFiles() was not declared static.
# o Made fromFiles() protected.
# 2008-03-11
# o BUG FIX: Calling readUnits(..., units=cdfUnits) twice, where 'cdfUnits'
#   was first created with, say, stratifyBy="pm" and then with
#   stratifyBy="pmmm", for the identical set of units, would give the same
#   results both times.  This was because the hash key list object for the
#   file cache where identical.
# o BUG FIX: getUnitIntensities() would not pass arguments '...', e.g.
#   'stratifyBy', to readUnits() for AffymetrixCdfFile.  Thanks Tim
#   Keighley, CSIRO, Sydney for reporting this.
# 2008-03-05
# o Updated message specifying that it looks for *sample* annotation files.
# 2008-02-28
# o Added getUnitGroupCellMap() to all AffymetrixCelSet classes.
# 2008-02-25
# o Now fromFiles() of AffymetrixCelSet uses the file cache, which speed up
#   the setup for extremely large data sets.  Note, currently the default
#   is 'force=TRUE' (for safety).
# 2008-01-30
# o Now as.character() for AffymetrixCelSet only reports time stamp if the
#   number of arrays in the data set is less that the threshold specified
#   in the 'aroma.affymetrix' settings.
# 2008-01-11
# o ROBUSTNESS: It was possible to set a non-compatible CDF when using
#   static fromFiles() of AffymetrixCelSet.
# o Added argument '.checkArgs' to setCdf() of AffymetrixCelSet, which now
#   validated the arguments once and not for every file.
# 2007-12-08
# o Added argument 'cdf' to static fromName() of AffymetrixCelSet.  When
#   using this argument, the 'chipType' argument is optional, and the
#   returned CEL set will be using the specified CDF.
# o Now fromFiles() of AffymetrixCelSet takes argument 'cdf' which can be
#   used to override the default CDF.  This way the default CDF must not
#   have to be in the search path.
# 2007-09-25
# o Now getAverageFile() will detect if an average file has been deleted
#   between calls and recalculate it.
# 2007-09-13
# o Added validate that the 'name' of an AffymetrixCelSet and
#   AffymetrixCelFile is of at least length one.
# 2007-09-06
# o Now setCdf() throws an (informative) error message whenever one tries
#   to use an ASCII CDF file. This behavior can be changed by setting
#   rule 'rules$allowAsciiCdfs' in the package settings.
# 2007-08-01
# o Renamed static fromName() of AffymetrixCelSet to byName().
# 2007-04-06
# o BUG FIX: fromFiles() of AffymetrixCelSet would give error "Exception:
#   Pathname not found: annotationData/samples" if that directory was
#   missing.  Now it is instead created.
# 2007-04-03
# o Now fromFiles() verifies that the set CDF is compatible with the CEL
#   files, otherwise the CDF of the first CEL file is used.
# 2007-04-02
# o Now sample annotation files are searched for in annotationData/samples/.
# 2007-03-28
# o Added argument 'cache=TRUE' to getUnitIntensities() and readUnits().
# 2007-03-24
# o BUG FIX: clearCache() did not clear the .readUnitsCache field.
# 2007-03-16
# o BUG FIX: getAverageFile() of AffymetrixCelSet would average the wrong
#   set of cells if argument 'indices' was different from NULL.
# 2007-03-06
# o Now attributes are set from SAF files in fromFiles().
# 2007-02-22
# o Fixed the warning about "'tzone' attributes are inconsistent". See
#   code of as.character() for explanation.
# o Now fromFiles() accepts argument 'chipType' to override any chip type
#   specified in the CEL headers. This is useful in case different CEL files
#   refers to different chip types, which can be the case for mixed
#   generations of CEL files.  Also added a scan of chip types.
# 2007-02-14
# o Added test for correct directory structure to fromFiles().  This will
#   enforce users to use the correct structure so that for instance the
#   name of the data set is correctly inferred.
# o Added argument 'onDuplicates' to fromFiles() so it is possible to
#   exclude duplicated CEL files.
# 2007-02-06
# o Added static method findByName() and fromName().
# 2007-01-15
# o Added 'classes=class(this)' to all "digest" keys.
# 2007-01-07
# o BUG FIX: In KS's update of getAverageFile() to support averaging
#   over other fields than intensities, argument 'indices' was missing
#   in the readCel() call making the function fail when processed chunk
#   by chunk. /HB
# 2007-01-05
# o Removed getSampleNames().
# 2006-12-01
# o Now as.character() reports the range of CEL header timestamps.
# 2006-11-07
# o Now getAverageFile() uses rowMedians() of R.native if available,
#   otherwise a local version utilizing apply(). Same for rowMads().
# 2006-10-24
# o Added getAverageLog() and getAverageAsinh().
# o Added transforms and anti-transforms g() and h() to getAverageFile().
# o Changed the defaults from mean to median, and sd to mad for
#   getAverageFile().
# o Added Rdoc comments to getAverageFile().
# 2006-10-10
# o Renamed rma and gcrma to rmaSummary and gcrmaSummary, to avoid clash
#   with existing functions.
# o Added gcrma() wrapper function.
# o Added rma() wrapper function.
# o Fixed a bug in getData() - default for argument "fields" contained "xy",
#   which is not a valid field (x, y are separate).
# 2006-10-02
# o Added getData().  Now getIntensities() works again and is just a wrapper
#   to getData().
# 2006-09-18
# o Now references to all requested average files are cached so it can
#   return the same object instead of creating a new one each time.
# 2006-09-16
# o Added getSiblings() to easily get other data sets for the same
#   samples.
# 2006-09-14
# o Added a read-buffer cache to readUnits() and getUnitIntensities().
# 2006-08-27
# o Added getAverageFile().
# 2006-08-26
# o Now getName() of a CEL set is inferred from the pathname:
#     path/to/<name>/chip_files/<"chip type">/
# 2006-08-21
# o Now AffymetrixCelSet inherits from AffymetrixFileSet.
# 2006-08-11
# o Added clearCache() which also clears the cache of all data file object.
# 2006-05-16
# o Redefined "[" to extract arrays.
# 2006-04-13
# o Added Rdoc comments for all methods.
# 2006-04-09
# o Now the read map is loaded automatically when fromFiles() used.
# 2006-03-30
# o Updated to new aroma.apd.
# 2006-03-18
# o Added argument 'subset' to calcAvgCellSignals() & normalizeQuantile().
# 2006-03-15
# o Now nbrOfCells() returns the number of cells for the first file only.
# o Now the fromFiles(static, ...) creates an object of the same class as
#   the static object.
# 2006-03-04
# o Added mapping functions.
# o Added writeApd().
# 2006-03-03
# o Added lapply().
# 2006-03-02
# o Updated to deal with AffymetrixDataFile object instead of CEL directly.
# 2006-02-21
# o Letting readCelUnits() transform signals improves speed substantially.
# o Making use of new multi-array readCelUnits().
# 2006-02-20
# o Created.
############################################################################
