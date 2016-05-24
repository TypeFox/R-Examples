###########################################################################/**
# @RdocClass RspComment
#
# @title "The RspComment class"
#
# \description{
#  @classhierarchy
#
#  An RspComment is an @see "RspConstruct" that represents an RSP comment,
#  which are of format \code{<\%-- ... --\%>}, \code{<\%--- ... ---\%>}
#  and so on.  They can also be so called "empty" RSP comments of format
#  \code{<\%-\%>}, \code{<\%--\%>}, \code{<\%---\%>} and so on.
# }
#
# @synopsis
#
# \arguments{
#   \item{str}{A @character string.}
#   \item{...}{Not used.}
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
setConstructorS3("RspComment", function(str=character(), ...) {
  extend(RspConstruct(str), "RspComment");
})


#########################################################################/**
# @RdocMethod getComment
#
# @title "Gets the comment"
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
#*/#########################################################################
setMethodS3("getComment", "RspComment", function(comment, ...) {
  as.character(comment);
})


setMethodS3("asRspString", "RspComment", function(object, ...) {
  body <- unclass(object);
  suffixSpecs <- attr(object, "suffixSpecs");
  fmtstr <- "<%%%s%s%%>";
  s <- sprintf(fmtstr, body, suffixSpecs);
  RspString(s);
})

##############################################################################
# HISTORY:
# 2013-03-26
# o Added asRspString() for RspComment.
# 2013-02-11
# o Added Rdoc help.
# 2013-02-09
# o Created.
##############################################################################
