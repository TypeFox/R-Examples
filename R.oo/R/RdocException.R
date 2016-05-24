###########################################################################/**
# @RdocClass RdocException
#
# @title "RdocException are thrown by the Rdoc compiler"
#
# \description{
#  @classhierarchy
#
#  @get "title" when it fails to generate a Rd file from an Rdoc comment.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Any arguments accepted by @Exception}.
#   \item{source}{Object specifying the source where the Rdoc error occured.
#     This is commonly a filename @character string.}.
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# \seealso{
#   For detailed information about exceptions see @Exception.
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setConstructorS3("RdocException", function(..., source=NULL) {
  extend(Exception(...), "RdocException",
    .source = source
  )
})



###########################################################################/**
# @RdocMethod as.character
#
# \title{Gets a character string representing of the RdocException}
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
# \keyword{error}
#*/###########################################################################
setMethodS3("as.character", "RdocException", function(x, ...) {
  # To please R CMD check
  this <- x;

  paste("[", getWhen(this), "] ", class(this)[1L], " in ", getSource(this),
                                          ": ", getMessage(this), sep = "");
})



###########################################################################/**
# @RdocMethod getSource
#
# \title{Gets the source of the exception}
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
#  Returns the source.
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
# \keyword{error}
#*/###########################################################################
setMethodS3("getSource", "RdocException", function(x, ...) {
  x$.source
})


############################################################################
# HISTORY:
# 2012-12-28
# o Replaced all data.class(obj) with class(obj)[1].
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2004-10-17
# o Added more Rdoc comments.
# 2003-04-28
# o Added the field source to refer to the source file in which the error
#   was found.
# 2003-04-12
# o Created.
############################################################################
