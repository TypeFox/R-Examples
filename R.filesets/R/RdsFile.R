###########################################################################/**
# @RdocClass RdsFile
# @alias loadObject.RdsFile
#
# @title "The RdsFile class"
#
# \description{
#  @classhierarchy
#
#  An RdsFile represents a binary file containing an R object
#  saved using the @see "base::saveRDS" function.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \seealso{
#   An object of this class is typically part of an @see "RdsFileSet".
# }
#*/###########################################################################
setConstructorS3("RdsFile", function(...) {
  extend(GenericDataFile(...), "RdsFile");
})


setMethodS3("loadObject", "RdsFile", function(this, ...) {
  loadRDS(this, ...);
})


###########################################################################/**
# @RdocGeneric loadRDS
# @alias loadRDS.default
# @alias loadRDS.RdsFile
#
# @title "Reads data from a RDS file"
#
# \description{
#  @get "title".
# }
#
# \usage{
#  @usage loadRDS,default
#  @usage loadRDS,RdsFile
# }
#
# \arguments{
#   \item{file}{A @character string, a @connection, or an @see "RdsFile"
#    specifying a RDS file/connection to be read.}
#   \item{...}{Additional arguments passed to @see "base::readRDS".}
# }
#
# \value{
#  Returns an R object.
# }
#
# @author
#
# \seealso{
#   @see "base::readRDS".
# }
#*/###########################################################################
setMethodS3("loadRDS", "default", function(file, ...) {
  readRDS(file, ...);
})

setMethodS3("loadRDS", "RdsFile", function(file, ...) {
  pathname <- getPathname(file);
  loadRDS(pathname, ...);
})


############################################################################
# HISTORY:
# 2013-11-27
# o Added a default loadRDS(), which is a name alias for base::readRDS().
#   Also added loadRDS() for RdsFile.
# 2013-11-20
# o Created.
############################################################################
