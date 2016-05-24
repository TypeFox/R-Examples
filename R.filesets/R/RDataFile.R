###########################################################################/**
# @RdocClass RDataFile
# @alias loadObject.RDataFile
# @alias loadObject
#
# @title "The RDataFile class"
#
# \description{
#  @classhierarchy
#
#  An RDataFile represents a binary file containing R objects
#  saved using the @see "base::save" function.
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
#   An object of this class is typically part of an @see "RDataFileSet".
# }
#*/###########################################################################
setConstructorS3("RDataFile", function(...) {
  extend(GenericDataFile(...), "RDataFile")
})


setMethodS3("loadObject", "RDataFile", function(this, drop=TRUE, ...) {
  env <- loadToEnv(this, ...)

  # Drop?
  if (drop) {
    names <- ls(envir=env, all.names=TRUE)
    nvars <- length(names)
    if (nvars == 0L) {
      return(NULL)
    } else if (nvars == 1L) {
      return(env[[names]])
    }
  }

  as.list(env)
})


###########################################################################/**
# @RdocGeneric loadToEnv
# @alias loadToEnv.RDataFile
#
# @title "Reads data from a RDS file"
#
# \description{
#  @get "title".
# }
#
# \usage{
#  @usage loadToEnv,RDataFile
# }
#
# \arguments{
#   \item{file}{A @character string, a @connection, or an @see "RDataFile"
#    specifying an RData file to be read.}
#   \item{...}{Additional arguments passed to @see "R.utils::loadToEnv".}
# }
#
# \value{
#  Returns an @environment.
# }
#
# @author
#
# \seealso{
#   @see "R.utils::loadToEnv".
# }
#*/###########################################################################
setMethodS3("loadToEnv", "RDataFile", function(file, ...) {
  pathname <- getPathname(file)
  loadToEnv(pathname, ...)
})


############################################################################
# HISTORY:
# 2015-01-12
# o Created from RdsFile.R.
############################################################################
