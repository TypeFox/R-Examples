###########################################################################/**
# @RdocDefault sourceAllRsp
#
# @title "Processes one or several RSP files"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pattern}{A filename pattern.}
#   \item{path}{The pathname of the directory to be search for RSP files.}
#   \item{extension}{The filename extension to be used for the output files.}
#   \item{outputPath}{The pathname of the directory where to save the 
#         output files.}
#   \item{overwrite}{If @FALSE, an error is thrown if an output file already
#     exists, otherwise not.}
#   \item{...}{Additional arguments passed to @see "sourceRsp".}
#   \item{envir}{An @environment to be the working environment of the 
#     servlets, i.e. where RSP variables and objects are stored.}
# }
#
# \value{
#   Returns (invisibly) a @character @list of pathnames of all processed
#   RSP files. 
# }
#
# \section{Exceptions}{
#   If an exception occurs while processing a file, the methods skips to
#   the next one and records the error.
# }
#
# @examples "../incl/sourceAllRsp.Rex"
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
setMethodS3("sourceAllRsp", "default", function(pattern="[.]rsp$", path=".", extension="html", outputPath=extension, overwrite=FALSE, ..., envir=parent.frame()) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pattern'
  pattern <- Arguments$getRegularExpression(pattern);

  # Argument 'path'
  path <- Arguments$getReadablePath(path=path, mustExist=TRUE);

  # Argument 'extension'
  extension <- Arguments$getCharacter(extension, length=1, nchar=c(1,64));

  # Argument 'outputPath'
  outputPath <- Arguments$getWritablePath(path=outputPath);

  # Argument 'envir'
  if (!is.environment(envir))
    throw("Argument 'envir' is not an environment: ", class(envir)[1]);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # MAIN
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  files <- list.files(pattern=pattern, path=path, all.files=TRUE, full.names=TRUE);
  res <- list();
  for (file in files) {
    # Replace current filename extension with the output extension
    output <- gsub("[.](.*)$", paste(".", extension, sep=""), basename(file));
    output <- filePath(outputPath, output);

    # Create a response object for the output file
    response <- FileRspResponse(file=output, overwrite=overwrite);

    # Process the RSP file
    tryCatch({
      sourceRsp(file=file, path=path, response=response, ..., envir=envir);
    }, error = function(ex) {
      cat("An exception occured while processing RSP file '", file, 
                                              "' (in '", path, "'):", sep="");
      print(ex);
      attr(file, "error") <<- ex;
    })
    res <- c(res, file);
  }

  invisible(res);
})

##############################################################################
# HISTORY:
# 2005-10-31
# o Added argument 'overwrite' to sourceAllRsp().
# 2005-08-02
# o Now the method continues to the next RSP file if an error is detected.
#   Errors are returned as an attribute of the file in the result list.
# 2005-08-01
# o Added Rdoc comments.
# o Added more argument validation.
# 2005-07-31
# o Created.
##############################################################################
