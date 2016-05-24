###########################################################################/**
# @RdocClass AffymetrixCelSetTuple
#
# @title "The AffymetrixCelSetTuple class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#      @see "aroma.core::AromaMicroarrayDataSetTuple".}
#   \item{.setClass}{The name of the class of the input set.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
#*/###########################################################################
setConstructorS3("AffymetrixCelSetTuple", function(..., .setClass="AffymetrixCelSet") {
  extend(AromaMicroarrayDataSetTuple(..., .setClass=.setClass),
                                                     "AffymetrixCelSetTuple");
})




setMethodS3("byPath", "AffymetrixCelSetTuple", function(static, path, ..., verbose=FALSE) {
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

  # Get the corresponding CEL set class
  className <- class(static)[1];
  className <- gsub("Tuple$", "", className);
  clazz <- Class$forName(className);

  verbose && enter(verbose, "Creating ", class(static)[1]);
  verbose && cat(verbose, "Path: ", path);

  # Get all subdirectories
  verbose && enter(verbose, "Scanning path for subdirectories");
  dirs <- list.files(path=path, full.names=TRUE);
  dirs <- dirs[sapply(dirs, FUN=isDirectory)];
  verbose && print(verbose, dirs);
  verbose && exit(verbose);

  # Check which corresponds to chip types, i.e. has CDF files
  verbose && enter(verbose, "Keeping those with names of chip types");
  dirs <- sapply(dirs, FUN=function(dir) {
    chipType <- basename(dir);
    pathname <- AffymetrixCdfFile$findByChipType(chipType);
    if (is.null(pathname))
      dir <- NA;
    dir;
  })
  dirs <- dirs[!is.na(dirs)];
  verbose && print(verbose, dirs);
  verbose && exit(verbose);

  # Define CEL sets for each directory
  verbose && enter(verbose, "Defining list of ", className, ":s");
  csList <- lapply(dirs, FUN=function(dir) {
    clazz$byPath(dir, ..., verbose=less(verbose));
  })
  names(csList) <- basename(names(csList));
  verbose && str(verbose, dirs);
  verbose && exit(verbose);

  verbose && exit(verbose);

  newInstance(static, csList, ...);
}, static=TRUE)



setMethodS3("getListOfCdfs", "AffymetrixCelSetTuple", function(this, ...) {
  csList <- getSets(this);
  lapply(csList, FUN=getCdf);
}, private=TRUE)


setMethodS3("getListOfUnitNamesFiles", "AffymetrixCelSetTuple", function(this, ...) {
  csList <- getSets(this);
  lapply(csList, FUN=getCdf);
}, private=TRUE)


setMethodS3("getListOfUnitTypesFiles", "AffymetrixCelSetTuple", function(this, ...) {
  csList <- getSets(this);
  lapply(csList, FUN=getCdf);
}, private=TRUE)




##############################################################################
# HISTORY:
# 2009-12-31
# o Replaced all getListOfSets() with getSets().
# 2009-07-08
# o Added getListOfUnitTypesFiles() to AffymetrixCelSetTuple.
# 2009-01-26
# o Added getListOfUnitNamesFiles() to AffymetrixCelSetTuple.
# o Removed obsolete clearCache() for AffymetrixCelSetTuple.
# 2008-07-30
# o BUG FIX: byPath() and getListOfCdfs() was defined for superclass
#   AromaMicroarrayDataSetTuple and not this class.
# 2008-05-16
# o Extract all platform-generic methods to AromaMicroarrayDataSetTuple,
#   which is in aroma.core.
# o Now getFullName() getTableOfArrays(), getArrays(), indexOfArrays(),
#   getArrayTuple(), and asMatrixOfFiles() passes down '...' (for the
#   purpose passing down argument 'translate').
# o Now getFullNames() passes down 'translate'.
# 2008-03-29
# o getTableOfArrays() of AromaMicroarrayDataSetTuple returned the incorrect
#   array indices for the 2nd chip type if different arrays in the two sets.
# 2008-03-11
# o Renamed getTuple() to getArrayTuple().
# 2007-03-29
# o Added asMatrixOfFiles().
# 2007-03-20
# o Now getArrays() returns a named list where the names are the result from
#   getFullNames().
# 2007-03-19
# o TODO: Handle replicated sample names. It is not clear how this should be
#   done.
# o Created from GladModel.R.
##############################################################################
