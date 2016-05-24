###########################################################################/**
# @RdocClass ArrayExplorer
#
# @title "The ArrayExplorer class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{csTuple}{An @see "AffymetrixCelSet" object.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
#*/###########################################################################
setConstructorS3("ArrayExplorer", function(csTuple=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csTuple':
  if (!is.null(csTuple)) {
    if (!inherits(csTuple, "AffymetrixCelSetTuple")) {
      csTuple <- AffymetrixCelSetTuple(csTuple);
    }
  }

  extend(Explorer(...), "ArrayExplorer",
    .csTuple = csTuple,
    .reporters = NULL
  )
})


setMethodS3("as.character", "ArrayExplorer", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  s <- c(s, sprintf("Tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Number of chip types: %d", nbrOfChipTypes(this)));
  s <- c(s, paste("Number of arrays:", nbrOfArrays(this)));
  colorMaps <- getColorMaps(this);
  if (length(colorMaps) == 0) {
    colorMaps <- "<not specified>";
  } else {
    colorMaps <- paste(colorMaps, collapse="; ")
  }
  s <- c(s, paste("Color maps:", colorMaps));
  s <- c(s, sprintf("Main path: %s", getMainPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)



setMethodS3("getListOfReporters", "ArrayExplorer", function(this, ...) {
  reporters <- this$.reporters;

  if (is.null(reporters)) {
    # No need to add the tags, because they are now automatically inferred from
    # the input set.  This is done by the new AffymetrixFileSetReporter superclass
    # of SpatialReporter.  /HB 2008-03-29
#    tags <- getTags(this);
    setTuple <- getSetTuple(this);
    csList <- getSets(setTuple);
    reporters <- lapply(csList, FUN=function(cs) {
      reporter <- SpatialReporter(cs, tags="*");
      reporter;
    });
    this$.reporters <- reporters;
  }

  reporters;
}, protected=TRUE);


setMethodS3("setAlias", "ArrayExplorer", function(this, ...) {
  NextMethod("setAlias");
  reporters <- getListOfReporters(this);
  lapply(reporters, FUN=setAlias, ...);
  invisible(this);
}, protected=TRUE)


setMethodS3("getAlias", "ArrayExplorer", function(this, ...) {
  reporters <- getListOfReporters(this);
  getAlias(reporters[[1]], ...);
}, protected=TRUE)


setMethodS3("getAsteriskTags", "ArrayExplorer", function(this, ...) {
  "";
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getDataSet
#
# @title "Gets the data set"
#
# \description{
#  @get "title" for which the explorer is displaying it results.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "AffymetrixCelSet".
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getDataSet", "ArrayExplorer", function(this, ...) {
  csTuple <- getSetTuple(this);
  csList <- getSets(csTuple);
  csList[[1]];  # AD HOC for now. /HB 2007-03-19
})

setMethodS3("getSetTuple", "ArrayExplorer", function(this, ...) {
  this$.csTuple;
})

setMethodS3("getNameOfInput", "ArrayExplorer", function(this, ...) {
  st <- getSetTuple(this);
  getName(st, ...);
}, protected=TRUE)


setMethodS3("getTagsOfInput", "ArrayExplorer", function(this, ...) {
  st <- getSetTuple(this);
  getTags(st, ...);
}, protected=TRUE)


setMethodS3("nbrOfChipTypes", "ArrayExplorer", function(this, ...) {
  nbrOfChipTypes(getSetTuple(this));
})

setMethodS3("getListOfUnitNamesFiles", "ArrayExplorer", function(this, ...) {
  getListOfUnitNamesFiles(getSetTuple(this));
}, private=TRUE)

setMethodS3("getListOfUnitTypesFiles", "ArrayExplorer", function(this, ...) {
  getListOfUnitTypesFiles(getSetTuple(this));
}, private=TRUE)


setMethodS3("getArraysOfInput", "ArrayExplorer", function(this, ...) {
  setTuple <- getSetTuple(this);
  getNames(setTuple);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod setArrays
#
# @title "Sets the arrays"
#
# \description{
#  @get "title" to be processed by the explorer.
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @character (or @integer) @vector of arrays to be
#      considered. If @NULL, all arrays of the data set are considered.}
#   \item{...}{Not used.}
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
setMethodS3("setArrays", "ArrayExplorer", function(this, ...) {
  tuple <- getSetTuple(this);
  idxs <- indexOf(tuple, ...);

  # Get the names of the arrays
  arrayNames <- getNames(tuple)[idxs];

  # Sanity check
  stopifnot(!any(duplicated(arrayNames)));

  this$.arrays <- arrayNames;
})



setMethodS3("addColorMap", "ArrayExplorer", function(this, ...) {
  reporters <- getListOfReporters(this);
  res <- lapply(reporters, FUN=addColorMap, ...);
  invisible(res);
})


setMethodS3("setColorMaps", "ArrayExplorer", function(this, ...) {
  reporters <- getListOfReporters(this);
  res <- lapply(reporters, FUN=setColorMaps, ...);
  invisible(res);
})

setMethodS3("getColorMaps", "ArrayExplorer", function(this, parsed=FALSE, ...) {
  # All reporters should have the same color maps
  reporter <- getListOfReporters(this)[[1]];
  getColorMaps(reporter, ...);
})


setMethodS3("updateSetupExplorerFile", "ArrayExplorer", function(this, ...) {
  setTuple <- getSetTuple(this);
  env <- new.env();
  env$reporters <- getListOfReporters(this);
  env$chipTypes <- getChipTypes(setTuple, fullname=TRUE);
  env$arrays <- getFullNames(setTuple);

  NextMethod("updateSetupExplorerFile", data=env);
}, private=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Generates image files, scripts and dynamic pages for the explorer"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{An optional @vector of arrays to be processed.
#      If @NULL, all arrays are considered.}
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns nothing.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "ArrayExplorer", function(this, arrays=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays)) {
  } else {
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays(this));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Generating ", class(this)[1], " report");

  # Setup HTML, CSS, Javascript files first
  setup(this, ..., verbose=less(verbose));

  # Get the array aliases
  aliases <- names(getArrays(this));
  if (!is.null(arrays)) {
    aliases <- aliases[arrays];
  }

  # Generate bitmap images
  reporters <- getListOfReporters(this);
  res <- lapply(reporters, FUN=function(reporter) {
    writeImages(reporter, arrays=arrays, aliases=aliases, ..., verbose=less(verbose));
  });

  # Update Javascript files
  updateSetupExplorerFile(this, ..., verbose=less(verbose));

  verbose && exit(verbose);
})


##############################################################################
# HISTORY:
# 2013-06-01
# o SPEEDUP: process() of ArrayExplorer no longer runs the garbage collector
#   at the end.
# 2012-10-18
# o CLEANUP: Drop all usage of 'version' in ArrayExplorer.
# 2012-03-06
# o Dropped setup() for ArrayExplorer, because now Explorer has one.
# o Replaced updateOnChipTypeJS() and updateOnLoadJS() with
#   updateSetupExplorerFile().
# 2009-07-08
# o Added getListOfUnitTypesFiles() for ArrayExplorer.
# 2009-05-19
# o Now testing for file permissions for creating/writing/updating files/dirs.
# 2009-01-26
# o Removed getListOfCdfs() from ArrayExplorer.
# 2008-08-19
# o Added argument 'arrays' to process() so that it is possible to specify
#   for which arrays images should be generated.
# 2008-06-03
# o BUG FIX: updateOnLoadJS() of ArrayExplorer did not use the fullname
#   chip type, cause an error in ArrayExplorer:s for tagged chip types.
#   Thanks to Maria Traka at BBSRC in UK for spotting this.
# 2008-03-29
# o BUG FIX: The ArrayExplorer would generate image files to a directory under
#   reports/<dataSet>/<tags>,<tags>/..., i.e. the tags where replicated.  This
#   is a bug introduced in the latest release.
#   Details: No need to add the tags, because they are now automatically
#   inferred from the input set.  This is done by the new
#   AffymetrixFileSetReporter superclass of SpatialReporter.
# 2007-08-09
# o Renamed updateSampleFile() to updateOnLoadJS().
# o Added updateOnChipTypeJS().
# 2007-03-20
# o Removed argument arrays from process().
# o Added setAlias() which also sets the alias on the reporters.
# o Added getAlias() to inherit the alias from the reporters.
# 2007-03-19
# o Class can now handle multiple chip types.
# o Class is now making use of the SpatialReporter class.
# o Now ChromosomeExplorer extends Explorer.
# o BUG FIX: setArrays() called indexOfArrays() instead of indexOf().
# 2007-02-28
# o BUG FIX: setColorMaps() gave "Error in addColorMap.ArrayExplorer(this,
#   colorMap, ...) : object "nbrOfColors" not found".
# 2007-02-08
# o Created from ChromosomeExplorer.R.
##############################################################################
