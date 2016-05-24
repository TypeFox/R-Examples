###########################################################################/**
# @RdocClass RspResponse
#
# @title "The RspResponse class"
#
# \description{
#  @classhierarchy
#
#  An abstract class that provides basic methods to write and flush output to
#  the generated document.
# }
#
# @synopsis
#
# \arguments{
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
setConstructorS3("RspResponse", function(...) {
  extend(Object(), "RspResponse",
    ...
  )
})




#########################################################################/**
# @RdocMethod write
# @alias write.FileRspResponse
# @alias write.HttpDaemonRspResponse
#
# @title "Writes an RSP response to the predefined output"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Objects to be pasted together and outputted.}
#   \item{collapse}{A @character string to be used to collapse the objects.}
#   \item{sep}{A @character string to separate the objects.}
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
#
# @keyword IO
#*/#########################################################################
setMethodS3("write", "RspResponse", abstract=TRUE);




#########################################################################/**
# @RdocMethod flush
# @alias flush.FileRspResponse
# @alias flush.HttpDaemonRspResponse
#
# @title "Flushes the response buffer"
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
#  Returns nothing.
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
setMethodS3("flush", "RspResponse", appendVarArgs=FALSE, abstract=TRUE)



###########################################################################/**
# @RdocMethod import
#
# @title "Imports the output from another RSP file"
#
# \description{
#  @get "title".
#  This is an internal methods called when processing the
#  \code{<\%@ import file="url"\%>} tag.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "sourceRsp".}
# }
#
# \value{
#   Writes the output from @see "sourceRsp" to the RSP response output file.
#   If an error occurs, the error message is written too.
# }
#
# @author
#
# \seealso{
#   @see "sourceRsp".
#   @seeclass
# }
#
# @keyword file
# @keyword IO
#*/###########################################################################
setMethodS3("import", "RspResponse", function(response, ...) {
  tryCatch({
    sourceRsp(..., response=response);
  }, error = function(ex) {
    error <- as.character(ex);
    write(response, error);
  })
})



##############################################################################
# HISTORY:
# 2006-07-04
# o Renamed from Response to RspResponse.
# 2005-11-30
# o Created from (File)RspResponse.  This is to be the new superclass.
##############################################################################
