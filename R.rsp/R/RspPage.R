###########################################################################/**
# @RdocClass RspPage
#
# @title "The RspPage class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{pathname}{A @character string.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
# @keyword internal
#*/###########################################################################
setConstructorS3("RspPage", function(pathname=NULL, ...) {
  # Argument 'pathname':
  pathname <- Arguments$getCharacter(pathname);

  extend(Object(), "RspPage",
    pathname = pathname,
    ...
  );
})



#########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of the directory of the current RSP file"
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
# @keyword IO
#*/#########################################################################
setMethodS3("getPath", "RspPage", function(this, ...) {
  getParent(this$pathname);
}, createGeneric=FALSE)



#########################################################################/**
# @RdocMethod getName
#
# @title "Gets the (base)name of the current RSP file"
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
# @keyword IO
#*/#########################################################################
setMethodS3("getName", "RspPage", function(this, ...) {
  basename(this$pathname);
}, createGeneric=FALSE)



#########################################################################/**
# @RdocMethod getAbsolutePath
#
# @title "Gets the absolute pathname to the current RSP file"
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
# @keyword IO
#*/#########################################################################
setMethodS3("getAbsolutePath", "RspPage", function(this, ...) {
  getAbsolutePath(this$pathname);
}, createGeneric=FALSE)



##############################################################################
# HISTORY:
# 2005-08-01
# o Added Rdoc comments.
# o Added getName() and getAbsolutePath().
# 2005-07-31
# o Created.
##############################################################################
