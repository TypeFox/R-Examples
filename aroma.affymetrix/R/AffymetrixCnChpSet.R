###########################################################################/**
# @RdocClass AffymetrixCnChpSet
#
# @title "The AffymetrixCnChpSet class"
#
# \description{
#  @classhierarchy
#
#  A AffymetrixCnChpSet object represents a set of AffymetrixCnChpFile:s
#  with \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "AffymetrixCnChpFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \seealso{
#   @see "AffymetrixCnChpFile".
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AffymetrixCnChpSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    reqFileClass <- "AffymetrixCnChpFile";
    lapply(files, FUN=function(df) {
      df <- Arguments$getInstanceOf(df, reqFileClass, .name="files");
    })
  } else if (inherits(files, "AffymetrixCnChpSet")) {
    return(as.AffymetrixCnChpSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  extend(AffymetrixFileSet(files=files, ...), "AffymetrixCnChpSet");
})



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the set"
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
setMethodS3("as.character", "AffymetrixCnChpSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  s <- c(s, sprintf("Tags: %s", tags));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  n <- length(this);
  s <- c(s, sprintf("Number of arrays: %d", n));
  names <- getNames(this);
  s <- c(s, sprintf("Names: %s [%d]", hpaste(names), n));
  s <- c(s, sprintf("Total file size: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)



setMethodS3("findByName", "AffymetrixCnChpSet", function(static, ..., paths="chpData(|,.*)/") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'paths':
  if (is.null(paths)) {
    paths <- eval(formals(findByName.AffymetrixCnChpSet)[["paths"]]);
  }

  NextMethod("findByName", paths=paths);
}, static=TRUE)


setMethodS3("byName", "AffymetrixCnChpSet", function(static, name, tags=NULL, chipType=NULL, cdf=NULL, paths=NULL, ...) {
  # Argument 'cdf':
  if (!is.null(cdf)) {
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");
  }

  # Argument 'chipType':
  if (is.null(chipType)) {
    if (!is.null(cdf)) {
      chipType <- getChipType(cdf, fullname=FALSE);  # Without tags
    } else {
      throw("Argument 'chipType' must be specified unless argument 'cdf' is specified.");
    }
  }

  suppressWarnings({
    path <- findByName(static, name, tags=tags, chipType=chipType, paths=paths, ...);
  })
  if (is.null(path)) {
    path <- file.path(paste(c(name, tags), collapse=","), chipType);
    throw("Cannot create ", class(static)[1], ".  No such directory: ", path);
  }

  suppressWarnings({
    byPath(static, path=path, cdf=cdf, ...);
  })
}, static=TRUE, protected=TRUE)


setMethodS3("byPath", "AffymetrixCnChpSet", function(static, path, pattern="[.](cnchp|CNCHP)$", cdf=NULL, checkChipType=is.null(cdf), ..., fileClass="AffymetrixCnChpFile", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Defining ", class(static)[1], " from files");

  # Call the "next" method
  # WORKAROUND: For unknown reasons it is not possible to specify
  # 'path=path' below, because it is already passed implicitly by
  # NextMethod() and if done, then argument 'recursive' to byPath() for
  # GenericDataFileSet will also be assign the value of 'path', e.g.
  # try byPath(AffymetrixCelSet(), "path/to/").  This seems to be related
  # to R-devel thread 'Do *not* pass '...' to NextMethod() - it'll do it
  # for you; missing documentation, a bug or just me?' on Oct 16, 2012.
  # [https://stat.ethz.ch/pipermail/r-devel/2012-October/065016.html]
  ##  set <- NextMethod("byPath", path=path, pattern=pattern, fileClass=fileClass, verbose=less(verbose));
  set <- NextMethod("byPath", pattern=pattern, fileClass=fileClass, verbose=less(verbose));

  verbose && cat(verbose, "Retrieved files: ", length(set));

  if (length(set) > 0) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Scan all CHP files for possible chip types
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Chip type according to the directory structure
    path <- getPath(set);
    chipType <- basename(path);
    verbose && cat(verbose,
                   "The chip type according to the path is: ", chipType);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Use the same CDF object for all CEL files.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for chip type '", chipType, "' inferred from path");
      cdf <- AffymetrixCdfFile$byChipType(chipType);
      verbose && exit(verbose);

      verbose && enter(verbose, "Check compatibility with 1st CEL file");
      verbose && cat(verbose, "Chip type: ", chipType);
      cf <- getOneFile(set);
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
  setCdf(set, cdf);
  verbose && exit(verbose);

  # Let the new CEL set update itself
  update2(set, verbose=less(verbose, 1));

  verbose && exit(verbose);

  set;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod as.AffymetrixCnChpSet
# @alias as.AffymetrixCnChpSet.list
# @alias as.AffymetrixCnChpSet.default
#
# @title "Coerce an object to an AffymetrixCnChpSet object"
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
#   Returns an @see "AffymetrixCnChpSet" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixCnChpSet", "AffymetrixCnChpSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixCnChpSet", "list", function(object, ...) {
  AffymetrixCnChpSet(object, ...);
})

setMethodS3("as.AffymetrixCnChpSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixCnChpSet object: ", mode(object));
})



setMethodS3("extractLogRatios", "AffymetrixCnChpSet", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract log ratios
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readMap <- NULL;
  data <- NULL;
  nbrOfArrays <- length(this);
  gcCount <- 0;
  for (kk in seq_len(nbrOfArrays)) {
    df <- this[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(df), nbrOfArrays));

    if (!is.null(readMap)) {
      setUnitReadMap(df, readMap=readMap);
    }

    dataKK <- extractLogRatios(df, units=units, ..., verbose=less(verbose, 5));

    if (is.null(readMap)) {
      readMap <- getUnitReadMap(df);
    }

    verbose && str(verbose, dataKK);
    if (is.null(data)) {
      dim <- c(length(dataKK), nbrOfArrays);
      dimnames <- list(NULL, getNames(this));
      naValue <- as.double(NA);
      data <- array(naValue, dim=dim, dimnames=dimnames);
    }
    data[,kk] <- dataKK;
    # Not needed anymore
    dataKK <- NULL;

    # Garbage collect?
    gcCount <- gcCount + 1;
    if (gcCount %% 10 == 0) {
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && exit(verbose);
  } # for (kk ...)

  # Drop singleton dimensions
  if (drop) {
    data <- drop(data);
  }

  verbose && cat(verbose, "Log ratios:");
  verbose && str(verbose, data);

  data;
})


setMethodS3("getCdf", "AffymetrixCnChpSet", function(this, ...) {
  getCdf(getOneFile(this), ...);
})


setMethodS3("setCdf", "AffymetrixCnChpSet", function(this, cdf, verbose=FALSE, ..., .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'cdf':
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Setting CDF for CEL set");
  verbose && print(verbose, cdf);

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


############################################################################
# HISTORY:
# 2011-02-24
# o Expanded the searched root paths to be chpData(|,.*)/
# 2009-08-12
# o Now findByName() of AffymetrixCnChpSet utilizes ditto of
#   AffymetrixCelSet, because its code was identical to the latter.
# 2008-08-22
# o Created.
############################################################################
