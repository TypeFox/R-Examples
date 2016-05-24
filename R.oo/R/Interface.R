###########################################################################/**
# @RdocClass Interface
#
# @title "The Interface class"
#
# \description{
#  @classhierarchy
#
#  This class represents a special set of classes whose purpose is to
#  provide methods (but not fields) shared by multiple different classes.
# }
#
# @synopsis
#
# \arguments{
#   \item{core}{The core value.}
#   \item{...}{Not used.}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
# @keyword internal
#*/###########################################################################
setConstructorS3("Interface", function(core=NA, ...) {
  this <- core;
  class(this) <- "Interface";
  this;
}, private=TRUE)


###########################################################################/**
# @RdocMethod extend
#
# @title "Extends another Interface class"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...className}{The name of new interface.}
#   \item{...}{Named values representing the fields of the new instance.}
# }
#
# \value{
#  Returns an Interface of class \code{className}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("extend", "Interface", function(this, ...className, ...) {
  class(this) <- unique(c(...className, class(this)));
  this;
})


###########################################################################/**
# @RdocMethod uses
# @alias uses
# @alias uses.character
#
# @title "Specifies that an object uses this Interface"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \value{
#  Returns a @character @vector of class names of Interface:s.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("uses", "Interface", function(this, ...) {
  res <- setdiff(class(this), "Interface");
  if (length(list(...)) > 0) {
    res <- c(list(res), list(uses(...)));

    # Order interfaces/classes
    names <- sort(unique(unlist(res, use.names=FALSE)));
    idxs <- integer(length(names));
    names(idxs) <- names;
    for (kk in seq_along(res)) {
      for (name in res[[kk]]) {
        idxs[name] <- kk;
      }
    }
    for (kk in seq_along(res)) {
      keep <- (idxs[res[[kk]]] == kk);
      res[[kk]] <- res[[kk]][keep];
    }
    res <- unlist(res, use.names=FALSE);
  }
  res;
})


setMethodS3("uses", "character", function(className, ...) {
  clazz <- Class$forName(className);
  obj <- newInstance(clazz);
  uses(obj, ...);
})



###########################################################################/**
# @RdocMethod as.character
#
# @title "Gets a character string representing the Interface"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("as.character", "Interface", function(x, ...) {
  # To please R CMD check
  this <- x;

  # Check if there are class "after" this one
  pos <- which("Interface" == class(this));
  isLast <- (pos == length(class(this)));
  if (isLast) {
    s <- paste(class(this), collapse=", ");
  } else {
    s <- NextMethod("as.character");
  }
  s;
}, private=TRUE)



###########################################################################/**
# @RdocMethod print
#
# @title "Prints an Interface"
#
# \description{
#  For all objects of class @see "Interface", this method will print the
#  value of \code{as.character()} of the object.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("print", "Interface", function(x, ...) {
  # To please R CMD check
  this <- x;

  print(as.character(this), ...);
})


###########################################################################/**
# @RdocMethod getFields
#
# @title "Returns NULL"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Ignored.}
# }
#
# \value{
#  Always returns @NULL.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("getFields", "Interface", function(...) { NULL }, private=TRUE)


############################################################################
# HISTORY:
# 2009-10-02
# o Added Rdoc comments.
# o Moved to R.oo.
# 2009-07-22
# o Now uses(...) takes multiple Interface classes.
# 2009-06-10
# o Added getFields() to Interface as an ad hoc solutions to avoid
#   print(<Interface>) throwing 'Error in UseMethod("getFields") : no
#   applicable method for "getFields"'.
# 2008-05-09
# o Added uses() for a character string.
# 2006-09-11
# o Added trial version of an Interface class.
############################################################################
