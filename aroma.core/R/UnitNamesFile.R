###########################################################################/**
# @RdocClass UnitNamesFile
#
# @title "The UnitNamesFile interface class"
#
# \description{
#  @classhierarchy
#
#  A UnitNamesFile provides methods for querying the unit names of
#  a given chip type.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "UnitAnnotationDataFile".}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("UnitNamesFile", function(...) {
  extend(UnitAnnotationDataFile(...), "UnitNamesFile");
})


setMethodS3("getUnitNames", "UnitNamesFile", abstract=TRUE);

setMethodS3("nbrOfUnits", "UnitNamesFile", function(this, ...) {
  length(getUnitNames(this));
})


###########################################################################/**
# @RdocMethod indexOf
#
# @title "Gets the indices of units by their names"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pattern}{A pattern to be used for identifying unit names of 
#      interest.  If @NULL, no regular expression matching is done.}
#   \item{names}{Names to be match exactly to the unit names.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @vector of @integers in [1,N] where N is the number of units
#  in the underlying annotation chip type file.
# }
#
# @author
#
# \seealso{
#   @seemethod "getUnitNames".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("indexOf", "UnitNamesFile", function(this, pattern=NULL, names=NULL, ...) {
  # Arguments 'pattern' & 'names':
  if (!is.null(pattern) && !is.null(names)) {
    throw("Only one of arguments 'pattern' and 'names' can be specified.");
  }

  if (is.null(pattern) && is.null(names)) {
    throw("Either argument 'names' or 'pattern' must be specified.");
  }

  if (!is.null(pattern)) {
    if (length(pattern) != 1) {
      throw("If specified, argument 'pattern' must be a single string. Did you mean to use argument 'names'?");
    }
    pattern <- Arguments$getRegularExpression(pattern);
    idxs <- grep(pattern, getUnitNames(this));
  } else if (!is.null(names)) {
    idxs <- match(names, getUnitNames(this));
  }

  idxs;
})


############################################################################
# HISTORY:
# 2009-07-08
# o CORRECTNESS: Added missing abstract getUnitNames().
# o CLEAN UP: Now UnitNamesFile inherits from UnitAnnotationDataFile.
# 2009-05-18
# o Now indexOf() for UnitNamesFile assert that exactly one of the 'pattern'
#   and 'names' arguments is given.  It also gives an informative error
#   message if 'pattern' is a vector.
# 2009-02-10
# o Added static byChipType() to UnitNamesFile.
# o Added a sanity check to getAromaUgpFile() of UnitNamesFile,
#   which asserts that the number of units in the located UGP file match
#   that of the data file. 
# 2009-01-26
# o Added getAromaUgpFile() to UnitNamesFile.
# 2008-07-21
# o Renamed UnitNamesInterface to UnitNamesFile.
# 2008-05-18
# o Created.
############################################################################
