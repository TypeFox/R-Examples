###########################################################################/**
# @RdocClass RspVariableDirective
# @alias RspStringDirective
# @alias RspNumericDirective
# @alias RspIntegerDirective
# @alias RspLogicalDirective
#
# @title "The RspVariableDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspVariableDirective is an @see "RspDirective" that causes the
#  RSP parser to assign the value of an attribute to an R object of
#  the same name as the attribute at the time of parsing.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Arguments passed to the constructor of @see "RspDirective".}
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
setConstructorS3("RspVariableDirective", function(value="variable", ...) {
  extend(RspDirective(value, ...), "RspVariableDirective")
})


setMethodS3("getInclude", "RspVariableDirective", function(object, ...) {
  attrs <- c("name", "content", "file", "default");
  has <- hasAttribute(object, attrs);
  names(has) <- attrs;
  if (!has["name"]) {
    return(FALSE);
  }
  if (!any(has[c("content", "file", "default")])) {
    return(TRUE);
  }
  FALSE;
})

setConstructorS3("RspStringDirective", function(value="string", ...) {
  extend(RspVariableDirective(value, ...), "RspStringDirective")
})

setConstructorS3("RspNumericDirective", function(value="numeric", ...) {
  extend(RspVariableDirective(value, ...), "RspNumericDirective")
})

setConstructorS3("RspIntegerDirective", function(value="integer", ...) {
  extend(RspNumericDirective(value, ...), "RspIntegerDirective")
})

setConstructorS3("RspLogicalDirective", function(value="logical", ...) {
  extend(RspVariableDirective(value, ...), "RspLogicalDirective")
})


###########################################################################/**
# @RdocClass RspMetaDirective
#
# @title "The RspMetaDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspMetaDirective is an @see "RspDirective" representing RSP metadata.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Arguments passed to the constructor of @see "RspStringDirective".}
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
setConstructorS3("RspMetaDirective", function(value="meta", ...) {
  extend(RspStringDirective(value, ...), "RspMetaDirective");
})



##############################################################################
# HISTORY:
# 2013-03-26
# o Now RspMetaDirective inherits RspStringDirective.
# 2013-03-15
# o Added requireAttributes() to RspDirective.
# o Added RSP meta directive.
# 2013-02-23
# o Added asRspString() for RspDirective and RspUnparsedDirective.
# 2013-02-22
# o Added RspUnparsedDirective.
# 2013-02-19
# o Added support for attribute 'text' of RspIncludeDirective:s.
# 2013-02-18
# o Added RspIfeqDirective, RspElseDirective, and RspEndifDirective.
# 2013-02-13
# o Added RspPageDirective.
# o Added 'language' attribute to RspEvalDirective.
# 2013-02-11
# o Added Rdoc help.
# 2013-02-09
# o Created.
##############################################################################
