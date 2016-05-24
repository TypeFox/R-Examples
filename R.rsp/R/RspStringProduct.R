###########################################################################/**
# @RdocClass RspStringProduct
#
# @title "The RspStringProduct class"
#
# \description{
#  @classhierarchy
#
#  An RspStringProduct is an @see RspProduct that represents an
#  RSP product in form of a @character string.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{@character strings passed to @see "RspProduct".}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setConstructorS3("RspStringProduct", function(...) {
  extend(RspProduct(...), "RspStringProduct");
})


#########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a plain character string representation of an RSP string product"
#
# \description{
#  @get "title".  All attributes including class have been dropped.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("as.character", "RspStringProduct", function(x, ...) {
  s <- unclass(x);
  attributes(s) <- NULL;
  s;
}, protected=TRUE)


setMethodS3("print", "RspStringProduct", function(x, ...) {
  print(as.character(x), ...);
}, protected=TRUE)




############################################################################
# HISTORY:
# 2013-02-13
# o Added RspProduct and RspFileProduct with corresponding
#   process() methods.
# o Created.
############################################################################
