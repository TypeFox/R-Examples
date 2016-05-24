###########################################################################/**
# @RdocDefault importRsp
#
# @title "Imports an RSP file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "sourceRsp".}
# }
#
# \value{
#   Returns the compile output of an RSP template as a @character string.
# }
#
# @author
#
# \seealso{
#   @see "sourceRsp".
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/########################################################################### 
setMethodS3("importRsp", "default", function(...) {
  input <- NULL; # Declare variable to please R CMD check R v2.6.0

  output <- textConnection("input", open="w", local=TRUE);
  on.exit(close(output));

  tryCatch({
    sourceRsp(..., response=FileRspResponse(file=output));
  }, error = function(ex) {
    error <- as.character(ex);
    input <<- paste(input, error, sep="");
    code <- ex$code;
    if (!is.null(code)) {
      code <- paste(code, collapse="\n", sep="");
      input <<- paste(input, code, sep="");
    }
  })

  input;
})


##############################################################################
# HISTORY:
# 2005-09-22
# o BUG FIX: sourceRsp() is no longer using argument 'output', but 'response'.
# 2005-09-18
# o Added Rdoc comments.
# 2005-07-31
# o Created.
##############################################################################
