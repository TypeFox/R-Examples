###########################################################################/**
# @RdocClass TransformReport
#
# @title "The TransformReport class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{inSet}{The input data set as an @see "AffymetrixCelSet".}
#   \item{outSet}{The output data set as an @see "AffymetrixCelSet".}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("TransformReport", function(inSet=NULL, outSet=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'inSet':
  if (!is.null(inSet)) {
    inSet <- Arguments$getInstanceOf(inSet, "AffymetrixCelSet");
    outSet <- Arguments$getInstanceOf(outSet, "AffymetrixCelSet");

    # Check for compatibility
#    if (!equals(getCdf(inSet), getCdf(outSet))) {
#      throw("Argument 'inSet' and 'outSet' have incompatible CDFs.");
#    }
  }

  extend(Object(), "TransformReport",
    .inSet = inSet,
    .outSet = outSet,
    .alias = NULL
  )
}, abstract=TRUE)



setMethodS3("getRootPath", "TransformReport", function(this, ...) {
  "reports";
}, private=TRUE)


setMethodS3("as.character", "TransformReport", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  ds <- getInputDataSet(this);
  s <- c(s, sprintf("Input data set: %s", getFullName(ds)));
  ds <- getOutputDataSet(this);
  s <- c(s, sprintf("Output data set: %s", getFullName(ds)));
  s <- c(s, sprintf("Number of arrays: %d (%.2fMB)",
                           length(ds), getFileSize(ds)/1024^2));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ds))));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output data set"
#
# \description{
#  @get "title", which is the same as the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "TransformReport", function(this, ...) {
  name <- getAlias(this);
  if (is.null(name)) {
    ds <- getOutputDataSet(this);
    name <- getName(ds);
  }
  name;
})

setMethodS3("getAlias", "TransformReport", function(this, ...) {
  this$.alias;
}, protected=TRUE)

setMethodS3("setAlias", "TransformReport", function(this, alias, ...) {
  if (!is.null(alias)) {
    alias <- Arguments$getCharacter(alias, nchar=c(1,Inf), length=c(1,1));
  }
  this$.alias <- alias;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which equals the tags of the input data set plus the tags
#  of this transformation.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "TransformReport", function(this, collapse=NULL, ...) {
  tags <- this$.tags;

  ds <- getOutputDataSet(this);
  tags <- getTags(ds, collapse=collapse);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output data set"
#
# \description{
#  @get "title", which is the name with comma separated tags.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "TransformReport", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of the output data set"
#
# \description{
#  @get "title".
#  If non-existing, then the directory is created.
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
#*/###########################################################################
setMethodS3("getPath", "TransformReport", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type
  ds <- getOutputDataSet(this);
  unf <- getUnitNamesFile(ds);
  chipType <- getChipType(unf, fullname=FALSE);

  # Image set
  set <- "transform";

  # The full path
  path <- filePath(rootPath, fullname, chipType, set);
  path <- Arguments$getWritablePath(path);

  path;
})



###########################################################################/**
# @RdocMethod getInputDataSet
#
# @title "Gets the source data set"
#
# \description{
#  @get "title" that is to be (or has been) transformed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet".
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getInputDataSet", "TransformReport", function(this, ...) {
  this$.inSet;
})



###########################################################################/**
# @RdocMethod getOutputDataSet
#
# @title "Gets the transformed data set"
#
# \description{
#  @get "title", if processed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet".
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getOutputDataSet", "TransformReport", function(this, ...) {
  this$.outSet;
})


setMethodS3("getCdf", "TransformReport", function(this, ...) {
  getCdf(getOutputDataSet(this));
})


setMethodS3("getUnitNamesFile", "TransformReport", function(this, ...) {
  dsOut <- getOutputDataSet(this);
  getUnitNamesFile(dsOut);
})

setMethodS3("getUnitTypesFile", "TransformReport", function(this, ...) {
  dsOut <- getOutputDataSet(this);
  getUnitTypesFile(dsOut);
})


setMethodS3("nbrOfArrays", "TransformReport", function(this, ...) {
  length(getOutputDataSet(this));
})

setMethodS3("seq", "TransformReport", function(this, ...) {
  seq_len(length(this));
})


setMethodS3("getYY", "TransformReport", function(this, array, transform=NULL, subset=1/8, field="intensities", ...) {
  inSet <- getInputDataSet(this);

  cdf <- getCdf(inSet);
  if (length(subset) == 1) {
    indices <- seq(from=1, to=nbrOfCells(cdf), length=subset*nbrOfCells(cdf));
    indices <- as.integer(indices);
  } else if (length(subset) > 1) {
    indices <- Arguments$getIndices(subset, max=nbrOfCells(cdf));
  } else {
    indices <- subset;
  }

  outSet <- getOutputDataSet(this);
  df1 <- inSet[[array]];
  df2 <- outSet[[array]];
  res <- list(
    array = array,
    y1 = getData(df1, indices=indices, ..., fields=field)[[field]],
    y2 = getData(df2, indices=indices, ..., fields=field)[[field]]
  );

  if (!is.null(transform)) {
    res$y1 <- transform(res$y1);
    res$y2 <- transform(res$y2);
  }

  res;
})



setMethodS3("plotXYCurve", "TransformReport", function(this, arrays=seq_along(this), lwd=2, col=arrays, xlim=c(0,65535), xlab=expression(y[1]), ylab=expression(y[2]), main=NULL, ..., add=FALSE, verbose=FALSE) {

  nbrOfArrays <- length(this);
  if (is.null(col)) {
    col <- seq_len(nbrOfArrays);
  } else {
    col <- rep(col, length.out=nbrOfArrays);
  }

  outSet <- getOutputDataSet(this);
  unf <- getUnitNamesFile(this);
  chipType <- getChipType(unf);
  for (kk in seq_along(arrays)) {
    array <- arrays[kk];

    df <- outSet[[array]];
    name <- getName(df);

    verbose && enter(verbose, sprintf("Array #%d ('%s')", kk, name));

    verbose && enter(verbose, "Retrieving data");
    suppressWarnings({
      yy <- getYY(this, array=array, ...);
    })
    verbose && str(verbose, yy);
    verbose && exit(verbose);

    if (is.null(main))
      main <- name;

    verbose && enter(verbose, "Plotting smooth (X,Y) curve");
    suppressWarnings({
      fit <- plotXYCurve(yy$y1, yy$y2, lwd=lwd, col=col[kk], xlim=xlim, xlab=xlab, ylab=ylab, ..., add=add);
    })
    verbose && exit(verbose);

    if (!add)
      stextChipType(chipType, line=-1);

    if (length(arrays) == 1)
      stextSize(df, size=length(yy$y1));

    add <- TRUE;

    # Not needed anymore
    yy <- NULL;

    # Garbage collect
    gc();

    verbose && exit(verbose);
  } # for (array in ...)

  box();

  invisible();
})


setMethodS3("plotXYCurveLog2", "TransformReport", function(this, xlim=c(0,16), xlab=expression(log[2](y[1])), ylab=expression(log[2]*(y[2])), ...) {
  plotXYCurve(this, transform=log2, xlim=xlim, xlab=xlab, ylab=ylab, ...);
})



setMethodS3("writeImages", "TransformReport", function(this, path=NULL, width=800, height=width, ..., skip=TRUE, verbose=FALSE) {
  pngDev <- findPngDevice(transparent=FALSE);

  unf <- getUnitNamesFile(this);
  chipType <- getChipType(unf, fullname=FALSE);
  rootPath <- getRootPath(this);
  name <- getName(this);
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  set <- "transform";
  path <- file.path(rootPath, name, tags, chipType, set);
  path <- Arguments$getWritablePath(path);

  outSet <- getOutputDataSet(this);
  nbrOfArrays <- length(outSet);

  verbose && enter(verbose, "Writing images for ", nbrOfArrays, " arrays");

  verbose && printf(verbose, "Image dimension: %.0fx%.0f\n", width, height);

  for (kk in seq_len(nbrOfArrays)) {
    df <- outSet[[kk]];
    fullname <- getFullName(df);
    verbose && enter(verbose, "Output CEL file: ", fullname);

    # Plot (log2(y),log2(y))
    tags <- c("YvY,log2");
    imgname <- paste(c(fullname, tags), collapse=",");
    filename <- sprintf("%s.png", imgname);
    pathname <- file.path(path, filename);

    verbose && cat(verbose, "Image pathname: ", pathname);
    if (!skip || !isFile(pathname)) {
      pngDev(pathname, width=width, height=height);
      tryCatch({
        plotXYCurveLog2(this, array=kk, dcol="#cccccc", ...);
      }, finally = {
        devDone();
      })
    }

    # Garbage collection
    if (kk %% 10 == 0)
      verbose && print(verbose, gc());

    verbose && exit(verbose);
  }

  # Garbage collection
  verbose && print(verbose, gc());

  verbose && exit(verbose);
})




setMethodS3("writeImageCombined", "TransformReport", function(this, path=NULL, width=800, height=width, ..., skip=TRUE, verbose=FALSE) {
  pngDev <- findPngDevice(transparent=FALSE);

  unf <- getUnitNamesFile(this);
  chipType <- getChipType(unf, fullname=FALSE);
  rootPath <- getRootPath(this);
  name <- getName(this);
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  set <- "transform";
  path <- file.path(rootPath, name, tags, chipType, set);
  path <- Arguments$getWritablePath(path);

  outSet <- getOutputDataSet(this);
  nbrOfArrays <- length(outSet);

  verbose && enter(verbose, "Writing combined image for ", nbrOfArrays, " arrays");

  verbose && printf(verbose, "Image dimension: %.0fx%.0f\n", width, height);

  # Plot (log2(y),log2(y))
  tags <- c("YvY,log2");
  imgname <- paste(c("all", tags), collapse=",");
  filename <- sprintf("%s.png", imgname);
  pathname <- file.path(path, filename);

  verbose && cat(verbose, "Image pathname: ", pathname);
  if (!skip || !isFile(pathname)) {
    pngDev(pathname, width=width, height=height);
    tryCatch({
      plotXYCurveLog2(this, ..., verbose=less(verbose));
    }, finally = {
      devDone();
    })
  }

  # Garbage collection
  verbose && print(verbose, gc());

  verbose && exit(verbose);
})




############################################################################
# HISTORY:
# 2011-11-07
# o Replaced dev.off() with devDone().
# 2009-07-08
# o Added getUnitTypesFile() for TransformReport.
# 2008-05-18
# o Made class less platform specific by utilizing UnitNamesFile interface.
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-02-06
# o Updated the path to <rootPath>/<dataSetName>/<tags>/<chipType>/<set>/.
# 2007-02-04
# o Created.
############################################################################
