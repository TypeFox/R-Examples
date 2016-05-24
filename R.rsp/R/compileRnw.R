###########################################################################/**
# @RdocDefault compileRnw
#
# @title "Compiles a Rnw file"
#
# \description{
#  @get "title".
#  The compiler used depends on the content type.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename, path}{The filename and (optional) path of the
#      document to be compiled.}
#   \item{...}{Additional arguments passed to the compiler function
#      used.}
#   \item{type}{A @character string specifying what content type of
#      Rnw file to compile.  The default (@NULL) is to infer the type
#      from the content of the file using @see "typeOfRnw".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname of the generated document.
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("compileRnw", "default", function(filename, path=NULL, ..., type=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- if (is.null(path)) filename else file.path(path, filename);
  if (!isUrl(pathname)) {
    pathname <- Arguments$getReadablePathname(pathname);
  }

  # Argument 'type':
  if (!is.null(type)) {
    type <- Arguments$getCharacter(type);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Compiling Rnw document");

  # Download URL?
  if (isUrl(pathname)) {
    verbose && enter(verbose, "Downloading URL");
    url <- pathname;
    verbose && cat(verbose, "URL: ", url);
    pathname <- downloadFile(url, verbose=less(verbose,50));
    verbose && cat(verbose, "Local file: ", pathname);
    verbose && exit(verbose);
  }

  # Now we can infer the type of Rnw, i.e. Sweave or knitr
  if (is.null(type)) {
    type <- typeOfRnw(pathname);
  }
  verbose && cat(verbose, "Type of Rnw file: ", type);

  if (type == "application/x-sweave") {
    pathnameR <- compileSweave(filename, path=path, ..., verbose=verbose);
  } else if (type == "application/x-knitr") {
    pathnameR <- compileKnitr(filename, path=path, ..., verbose=verbose);
  } else if (type == "application/x-asciidoc-noweb") {
    pathnameR <- compileAsciiDocNoweb(filename, path=path, ..., verbose=verbose);
  } else {
    throw("Unknown value of argument 'type': ", type);
  }

  verbose && exit(verbose);

  pathnameR;
}) # compileRnw()


############################################################################
# HISTORY:
# 2013-03-29
# o Added support for AsciiDoc Rnw:s.
# 2013-01-20
# o Created from compileSweave.R.
############################################################################
