#########################################################################/**
# @RdocDefault updateApdHeader
#
# @title "Updates the header of an Affymetrix probe data (APD) file"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The filename of the APD file.}
#   \item{path}{The path to the APD file.}
#   \item{...}{A set of named header values to be updated/added to the
#      header.  A value of @NULL will be removed from the current header.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns (invisibly) the pathname of the file updated.
# }
#
# @author
#
# \examples{\dontrun{#See ?createApd for an example.}}
#
# \seealso{
#   @see "createApd" and @see "updateApd".
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("updateApdHeader", "default", function(filename, path=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(filename, path=path,
                                                          mustExists=TRUE);

  # Argument 'nbrOfCells':
  nbrOfCells <- Arguments$getDouble(nbrOfCells, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get current header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  header <- readApdHeader(filename, path=path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create APD header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the header elements to be set
  args <- list(...);

  # Do not update reserved header elements
  reserved <- c("creator", "dataType", "bytesPerCell", "RStorageMode");
  excl <- (names(args) %in% reserved);
  args <- args[!excl];

  if (length(args) == 0) {
    return(invisible(pathname));
  }

  # Add to current header
  for (kk in seq(along=args)) {
    name <- names(args)[kk];
    value <- args[[kk]];
    header[[name]] <- value;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Wrap up the APD header in the file vector header comments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  comments <- c();
  for (kk in seq(along=header)) {
    key <- names(header)[kk];
    value <- header[[kk]];
    valueStr <- paste(key, "=", value, sep="");
    comments <- c(comments, valueStr);
  }
  comments <- paste(comments, collapse="\n");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the file vector
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  apd <- FileVector(pathname);
  on.exit(close(apd));

  setComments(apd, comments);

  invisible(pathname);
})


############################################################################
# HISTORY:
# 2006-03-14
# o Created.
############################################################################
