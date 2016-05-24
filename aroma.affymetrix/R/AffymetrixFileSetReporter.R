###########################################################################/**
# @RdocClass AffymetrixFileSetReporter
#
# @title "The AffymetrixFileSetReporter class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{set}{An @see "AffymetrixFileSet" object.}
#   \item{...}{Arguments passed to @see "GenericReporter".}
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
setConstructorS3("AffymetrixFileSetReporter", function(set=NULL, ..., .setClass="AffymetrixFileSet") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'set':
  if (!is.null(set)) {
    set <- Arguments$getInstanceOf(set, .setClass);
  }


  extend(GenericReporter(...), "AffymetrixFileSetReporter",
    .set = set
  )
})


setMethodS3("getFileSet", "AffymetrixFileSetReporter", function(this, ...) {
  this$.set;
}, protected=TRUE)


setMethodS3("getInputName", "AffymetrixFileSetReporter", function(this, ...) {
  set <- getFileSet(this);
  getName(set);
})

setMethodS3("getInputTags", "AffymetrixFileSetReporter", function(this, ...) {
  set <- getFileSet(this);
  getTags(set);
})



##############################################################################
# HISTORY:
# 2007-04-12
# o Now inheriting from the new GenericReporter class.
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-03-19
# o Created from Explorer.R.
##############################################################################
