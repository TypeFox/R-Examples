###########################################################################/**
# @RdocFunction compileRdoc
#
# @title "Compiles the Rdoc help comments in a package"
#
# \description{
#  @get "title" by parsing the source code files and writes Rd help files.
# }
#
# @synopsis
#
# \arguments{
#   \item{pkgname}{A @character string specifying the package name.
#     If @NULL, the package name is inferred from the DESCRIPTION file
#     that is automatically search for in the subdirectories.}
#   \item{path}{A @character string specifying the path of the
#     package source directory.}
#   \item{...}{Additional arguments passed to \code{Rdoc\$compile()}.}
#   \item{verbose}{If @TRUE, verbose output is printed, otherwise not.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \details{
#   To compile all Rdoc comments in a package 'PkgA' that is located
#   under the current working directory (e.g. \code{PkgA/R/*.R}) from
#   the system command line, do:
#   \preformatted{
#     Rscript -e R.oo::compileRdoc()
#   }
# }
#
# @author
#
# \seealso{
#   See the @see "Rdoc" class.
# }
#
# @keyword documentation
# @keyword IO
# @keyword internal
#*/###########################################################################
compileRdoc <- function(pkgname=NULL, path=pkgname, ..., verbose=TRUE) {
  require("R.oo") || stop("Package not loaded: R.oo");

  # Infer package name from DESCRIPTION?
  if (is.null(pkgname)) {
    dirs <- list.files();
    dirs <- dirs[file_test("-d", dirs)];
    pathnames <- file.path(dirs, "DESCRIPTION");
    pathnames <- pathnames[file_test("-f", pathnames)];
    if (length(pathnames) == 0L) {
      throw("Failed to infer package name, since no */DESCRIPTION file was found.");
    }
    if (length(pathnames) > 1L) {
      throw("Failed to infer unique package name, because more than one */DESCRIPTION file was found: ", paste(sQuote(pathnames), collapse=", "));
    }
    pathname <- pathnames[1L];
    pi <- read.dcf(file=pathname);
    pkgname <- pi[,"Package", drop=TRUE];
    if (length(pkgname) == 0L) {
      throw("Failed to infer package name.  No 'Package' was listed in ", sQuote(pathname), ".");
    }
    if (length(pkgname) > 1L) {
      throw("Failed to infer package name.  More than one 'Package' were listed in ", sQuote(pathname), ": ", paste(sQuote(pkgname), collapse=", "));
    }
  }

  if (is.null(path)) {
    path <- pkgname;
  }

  if (!file_test("-d", path)) {
    throw("No such package directory: ", path);
  }
  pathR <- file.path(path, "R")
  if (!file_test("-d", pathR)) {
    throw("No such package R/ directory: ", pathR);
  }


  require(pkgname, character.only=TRUE) || throw("Package not loaded: ", pkgname);

  opwd <- setwd(pathR);
  on.exit(setwd(opwd));

  Rdoc$compile(..., verbose=verbose);
} # compileRdoc()


############################################################################
# HISTORY:
# 2013-05-30
# o Added argument 'path' to compileRdoc().
# 2013-04-03
# o CLEANUP: compileRdoc() was outputting search() and sessionInfo().
# 2013-03-08
# o Created.
############################################################################
