###########################################################################/**
# @RdocClass DChipDcpSet
#
# @title "The DChipDcpSet class"
#
# \description{
#  @classhierarchy
#
#  A DChipDcpSet object represents a set of DChip DCP files
#  for \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "DChipDcpFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \seealso{
#   @see "DChipDcpFile".
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("DChipDcpSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    reqFileClass <- "DChipDcpFile";
    lapply(files, FUN=function(df) {
      df <- Arguments$getInstanceOf(df, reqFileClass, .name="files");
    })
  } else if (inherits(files, "DChipDcpSet")) {
    return(as.DChipDcpSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  extend(AffymetrixFileSet(files=files, ...), "DChipDcpSet");
})



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the DChip CHP set"
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
setMethodS3("as.character", "DChipDcpSet", function(x, ...) {
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



setMethodS3("findByName", "DChipDcpSet", function(static, ..., paths=c("rawData(|,.*)/", "probeData(|,.*)/")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'paths':
  if (is.null(paths)) {
    paths <- eval(formals(findByName.DChipDcpSet)[["paths"]]);
  }

  NextMethod("findByName", paths=paths);
}, static=TRUE, protected=TRUE)


setMethodS3("byName", "DChipDcpSet", function(static, name, tags=NULL, chipType, paths=NULL, ...) {
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

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



setMethodS3("byPath", "DChipDcpSet", function(static, path="rawData/", pattern="[.](dcp|DCP)$", ..., fileClass="DChipDcpFile", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  verbose && cat(verbose, "Retrieved files: ", length(this));


  if (length(this) > 0) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Scan all CHP files for possible chip types
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Chip type according to the directory structure
    path <- getPath(this);
    chipType <- basename(path);
    verbose && cat(verbose,
                   "The chip type according to the path is: ", chipType);
  }

  verbose && exit(verbose);

  this;
}, static=TRUE, protected=TRUE)




###########################################################################/**
# @RdocMethod as.DChipDcpSet
# @alias as.DChipDcpSet.list
# @alias as.DChipDcpSet.default
#
# @title "Coerce an object to an DChipDcpSet object"
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
#   Returns an @see "DChipDcpSet" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.DChipDcpSet", "DChipDcpSet", function(object, ...) {
  object;
})

setMethodS3("as.DChipDcpSet", "list", function(object, ...) {
  DChipDcpSet(object, ...);
})

setMethodS3("as.DChipDcpSet", "default", function(object, ...) {
  throw("Cannot coerce object to an DChipDcpSet object: ", mode(object));
})



setMethodS3("extractTheta", "DChipDcpSet", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
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
  # Extract the thetas
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- NULL;
  nbrOfArrays <- length(this);
  gcCount <- 0;
  for (kk in seq_len(nbrOfArrays)) {
    df <- this[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(df), nbrOfArrays));

    dataKK <- extractTheta(df, units=units, ..., verbose=less(verbose, 5));
    verbose && str(verbose, dataKK);
    if (is.null(data)) {
      dim <- c(nrow(dataKK), ncol(dataKK), nbrOfArrays);
      dimnames <- list(NULL, NULL, getNames(this));
      naValue <- as.double(NA);
      data <- array(naValue, dim=dim, dimnames=dimnames);
    }
    data[,,kk] <- dataKK;
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

  verbose && cat(verbose, "Thetas:");
  verbose && str(verbose, data);

  data;
})


############################################################################
# HISTORY:
# 2013-07-28
# o ROBUSTNESS: byPath() for DChipDcpSet was not decleared static.
# 2011-02-24
# o Expanded the searched root paths to be rawData(|,.*)/ and
#   probeData(|,.*)/.
# 2009-08-12
# o Now findByName() of DChipDcpSet utilizes ditto of AffymetrixCelSet,
#   because its code was identical to the latter.
# 2008-08-20
# o Added extractTheta().
# 2008-07-21
# o Now findByName() assert that the data set name is not empty.
# 2008-05-09
# o Now DChipDcpSet inherits from GenericDataFileSet.
# 2008-05-08
# o If paths=NULL in findByName(), it becomes the default argument value.
# 2008-01-30
# o Created.
############################################################################
