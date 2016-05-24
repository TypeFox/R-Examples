###########################################################################/**
# @RdocDefault readDataFrame
#
# @title "Reads data from a tabular file"
#
# \description{
#  @get "title" or a set of such files.
# }
# 
# @synopsis
#
# \arguments{
#   \item{filename, path}{@character @vector specifying one or more files to
#    be read.}
#   \item{...}{Additional arguments passed to either
#      (i) \code{\link[=readDataFrame.TabularTextFile]{readDataFrame}} 
#          for class @see "TabularTextFile", or
#      (ii) \code{\link[=readDataFrame.TabularTextFileSet]{readDataFrame}} 
#          for class @see "TabularTextFileSet",
#     depending on whether one or multiple files are read.
#   }
# }
#
# \value{
#  Returns a @data.frame.
# }
#
# \details{
#   When reading multiple files at once, first each file is read into
#   a @data.frame, and then these @data.frames are (by default) merged into
#   one @data.frame using @see "base::rbind".  This requires that the
#   same set of columns are read for each file.  Which columns to read
#   can be controlled by specifying their names in 
#   argument \code{colClasses}.  To change how the @data.frames are 
#   merged, use argument \code{combineBy}.
#   For more information, follow the help on the above to
#   \code{readDataFrame()} help links.
# }
#
# @examples "../incl/readDataFrame.Rex"
#
# @author
#
# \seealso{
#   @see "utils::read.table".
#   For further details, see classes @see "TabularTextFile" and 
#   @see "TabularTextFileSet".
# }
#*/###########################################################################
setMethodS3("readDataFrame", "default", function(filename, path=NULL, ...) {
  # Argument 'filename' and 'path':
  if (length(filename) == 0L) {
    throw("Argument 'filename' must not be empty.");
  }
  if (length(path) > 0L) {
    path <- rep(path, length.out=length(filename));
  }


  # Setup TabularTextFile or TabularTextFileSet
  if (length(filename) == 1L) {
    db <- TabularTextFile(filename, path=path, .verify=FALSE);
  } else if (length(filename) > 1L) {
    if (length(path) > 0L) {
      pathnames <- file.path(path, filename);
    } else {
      pathnames <- filename;
    }
    dfList <- lapply(pathnames, FUN=TabularTextFile, .verify=FALSE);
    db <- TabularTextFileSet(dfList);
  }

  readDataFrame(db, ...);
})



############################################################################
# HISTORY:
# 2013-01-17
# o Now readDataFrame() can read multiple files.
# 2013-01-16
# o Created.
############################################################################
