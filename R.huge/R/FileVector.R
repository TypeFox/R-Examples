###########################################################################/**
# @RdocClass FileVector
# @alias FileByteVector
# @alias FileShortVector
# @alias FileIntegerVector
# @alias FileFloatVector
# @alias FileDoubleVector
#
# @title "Class representing a persistent vector stored on file"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "AbstractFileArray".}
#  \item{length}{The number of elements in the vector.}
#  \item{names}{Optional element names.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
#
# }
#
# \details{
#   The purpose of this class is to be able to work with large vectors
#   in \R without being limited by the amount of memory available.
#   Data is kept on the file system and elements are read and written
#   whenever queried.
#
#   For more details, @see "AbstractFileArray".
# }
#
# \section{Supported data types}{
#   The following subclasses implement support for various data types:
#   \itemize{
#    \item \code{FileByteVector} (1 byte per element), 
#    \item \code{FileShortVector} (2 bytes per element), 
#    \item \code{FileIntegerVector} (4 bytes per element), 
#    \item \code{FileFloatVector} (4 bytes per element), and
#    \item \code{FileDoubleVector} (8 bytes per element).
#   }
# }
#
# @examples "../incl/FileVector.Rex"
#
# @author
#
# @visibility public
#*/########################################################################### 
setConstructorS3("FileVector", function(..., length=NULL, names=NULL) {
  dim <- length;
  dimnames <- list(names);

  extend(AbstractFileArray(..., dim=dim, dimnames=dimnames), "FileVector");
})



###########################################################################/**
# @RdocMethod names
#
# @title "Gets the element names of a file vector"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("names", "FileVector", function(x, ...) {
  # To please R CMD check.
  this <- x;

  dimnames(this)[[1]];
})




###########################################################################/**
# @RdocMethod "["
#
# @title "Gets a subset of elements of a file vector"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{i}{A @numeric @vector or element indices.}
# }
#
# \value{
#  Returns a @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "[<-".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("[", "FileVector", function(this, i=NULL) {
  if (is.null(i)) {
    values <- readAllValues(this);
  } else {
    values <- readValues(this, indices=i);
  }

  names <- names(this);
  if (length(names) > 0)
    names(values) <- names;

  values;
}) # "["()




###########################################################################/**
# @RdocMethod "[<-"
#
# @title "Assigns values to a subset of elements of a file vector"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{i}{A @numeric @vector or element indices.}
#   \item{value}{Values to be assigned to the selected elements.}
# }
#
# \value{
#  Returns itself.
# }
#
# @author
#
# \seealso{
#   @seemethod "[".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("[<-", "FileVector", function(this, i=NULL, value) {
  if (is.null(i)) {
    writeAllValues(this, values=value);
  } else {
    writeValues(this, indices=i, values=value);
  }

  this;
})



############################################################################
# HISTORY:
# 2006-05-09
# o Added Rdoc comments.
# 2006-02-27
# o Now inheriting from class AbstractFileArray.
# o Created from FileMatrix.R.
############################################################################ 
