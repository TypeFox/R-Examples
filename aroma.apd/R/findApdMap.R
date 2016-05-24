########################################################################/**
# @RdocDefault findApdMap
#
# @title "Search for APD map files in multiple directories"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{mapType}{A @character string of the map type to search for.}
#  \item{paths}{A @character @vector of paths to be searched.
#    The current directory is always searched at the beginning.
#    If @NULL, default paths are searched.  For more details, see below.}
#  \item{pattern}{A regular expression file name pattern to match.}
#  \item{...}{Additional arguments passed to @see "R.utils::findFiles".}
# }
#
# \value{
#  Returns a @vector of the full pathnames of the files found.
# }
#
# \details{
#   Note, the current directory is always searched at the beginning.
#   This provides an easy way to override other files in the search path.
#   If \code{paths} is @NULL, then a set of default paths are searched.
#   The default search path is consituted of:
#   \enumerate{
#    \item \code{"."}
#    \item \code{getOption("AFFX_APD_PATH")}
#    \item \code{Sys.getenv("AFFX_APD_PATH")}
#   }
#
#   One of the easiest ways to set system variables for \R is to
#   set them in an \code{.Renviron} file, see @see "base::Startup"
#   for more details.
# }
#
# @author
#
# @keyword file
# @keyword IO
#**/#######################################################################
setMethodS3("findApdMap", "default", function(mapType=NULL, paths=NULL, pattern="[.](a|A)(p|P)(m|M)$", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'paths':
  if (is.null(paths)) {
    paths <- paste(".",
                   getOption("AFFX_APD_PATH"),
                   Sys.getenv("AFFX_APD_PATH"),
             sep=";", collapse=";");
  }

  # Argument 'mapType':
  if (!is.null(mapType)) {
    pattern <- paste(mapType, pattern, sep="");
  }

  findFiles(pattern=pattern, paths=paths, ...);
})


############################################################################
# HISTORY:
# 2006-03-14
# o Created from findCdf.R now in the affxparser package.
############################################################################
