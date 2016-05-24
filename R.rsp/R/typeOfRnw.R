###########################################################################/**
# @RdocDefault typeOfRnw
#
# @title "Checks whether an Rnw file is a Sweave or a knitr file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename, path}{The filename and (optional) path of the Rnw file.}
#   \item{default}{A @character string specifying the default result.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   TBA.
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("typeOfRnw", "default", function(filename, path=NULL, default="application/x-sweave", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- if (is.null(path)) filename else file.path(path, filename);
  if (!isUrl(pathname)) {
    withoutGString({
      pathname <- Arguments$getReadablePathname(pathname);
    })
  }

  # Argument 'default':
  default <- Arguments$getCharacter(default);


  # Read content
  bfr <- readLines(pathname, warn=FALSE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for knitr-specific commands
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (any(regexpr("opts_chunk$set(", bfr, fixed=TRUE) != -1L)) {
    return("application/x-knitr");
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for Sweave-specific commands
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (any(regexpr("\\SweaveOpts(", bfr, fixed=TRUE) != -1L)) {
    return("application/x-sweave");
  }

  default;
}) # typeOfRnw()


############################################################################
# HISTORY:
# 2014-04-18
# o Now typeOfRnw() supports URLs.
# 2013-01-20
# o Created.
############################################################################
