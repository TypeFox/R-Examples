###########################################################################/**
# @RdocDefault rspToHtml
#
# @title "Compiles an RSP file to an HTML file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{file}{The filename of the RSP file to be compiled.}
#   \item{path}{An optional path to the RSP file.}
#   \item{outFile}{The filename of the output file.  
#     If @NULL, a default output file is used.}
#   \item{outPath}{An optional path to the output file.}
#   \item{extension}{The filename extension of the default output file.}
#   \item{overwrite}{If @TRUE, an existing output file is overwritten.}
#   \item{...}{Additional arguments passed to @see "sourceRsp".}
# }
#
# \value{
#   Returns the pathname to the generated document.
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
setMethodS3("rspToHtml", "default", function(file=NULL, path=NULL, outFile=NULL, outPath=NULL, extension="html", overwrite=TRUE, ...) {
  # Argument 'file' and 'path':
  pathname <- Arguments$getReadablePathname(file, path=path, mustExist=FALSE);

  # Argument 'extension':
  extension <- Arguments$getCharacter(extension);

  # Generate the name of the output document
  if (is.null(outFile)) {
    repl <- paste(".", extension, sep="");
    outFile <- gsub("[.](r|R)(s|S)(p|P)$", repl, basename(pathname));
  }

  # Generate the response object
  response <- FileRspResponse(file=outFile, path=outPath, overwrite=overwrite);

  # Assure that we do not overwrite the RSP file
  if (getAbsolutePath(pathname) == getAbsolutePath(response)) {
    throw("Cannot compile RSP document. Pathname of output document is identical to the pathname of the RSP file: ", pathname);
  }

  # Compile
  sourceRsp(file=pathname, response=response, ...);

  invisible(getAbsolutePath(response));
}, private=TRUE) # rspToHtml()


############################################################################
# HISTORY:
# 2007-01-07
# o BUG FIX: Argument 'path' was set to the directory of the 'file'.
# 2006-08-06
# o Created for conveniency.
############################################################################


