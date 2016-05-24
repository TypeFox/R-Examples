###########################################################################/**
# @RdocDefault extend
#
# @title "Extends a object"
#
# \description{
#   via a mechanism known as "parasitic inheritance".
#   Simply speaking this method "extends" the class of an object. What is actually
#   happening is that it creates an instance of class name \code{...className},
#   by taking another object and add \code{...className} to the class
#   list and also add all the named values in @... as attributes.
#
#   The method should be used by the constructor of a class and nowhere else.
# }
#
# @synopsis
#
# \arguments{
#   \item{this}{Object to be extended.}
#   \item{...className}{The name of new class.}
#   \item{...}{Attribute fields of the new class.}
# }
#
# \value{
#  Returns an object of class \code{...className}.
# }
#
# @author
#
# @examples "../incl/extend.default.Rex"
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("extend", "default", function(this, ...className, ...) {
  # Add class '...className' to the class vector
  class(this) <- unique(c(...className, class(this)));

  # Add attributes
  args <- list(...);
  names <- names(args);
  for (kk in seq_along(args)) {
    name <- names[kk];
    if (name == "class")
      throw("Trying to set class attribute: ", as.character(this));
    arg <- args[[kk]];
    attr(this, name) <- arg;
  }

  # Return modified object
  this;
}, conflict="quiet")


######################################################################
# HISTORY:
# 2006-02-09
# o Using 'conflict="quiet"' when creating extend.default() to avoid
#   the warning that the method already exists.  One less warning.
# 2005-06-14
# o This will generate a warning that "Method already existed and was
#   overwritten", but this is ok.
# 2005-06-08
# o Created.
######################################################################
