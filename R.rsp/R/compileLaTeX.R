###########################################################################/**
# @RdocDefault compileLaTeX
#
# @title "Compiles a LaTeX file"
#
# \description{
#  @get "title" to either PDF or DVI.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename, path}{The filename and (optional) path of the
#      LaTeX document to be compiled.
#      Only *.tex and *.ltx filename extensions are allowed.}
#   \item{format}{A @character string specifying the output format.}
#   \item{clean, quiet, texinputs}{Additional arguments passed to
#      @see "tools::texi2dvi".}
#   \item{...}{Not used.}
#   \item{outPath}{The output and working directory.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname of the generated (PDF or DVI) document.
# }
#
# \section{Supported filename extensions}{
#   Internally @see "tools::texi2dvi" is used, which in turn uses
#   \code{Sys.which("texi2dvi")} if available.  Most known implementation
#   of the latter will only recognize LaTeX documents with filename
#   extensions *.tex and *.ltx (case sensitive).  (Any other filenames
#   will be compiled with 'texinfo', which is not a LaTeX compiler.)
# }
#
# @author
#
# \seealso{
#   Internally, @see "tools::texi2dvi" is used.
#   To compile Sweave LaTeX documents, @see "compileSweave".
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("compileLaTeX", "default", function(filename, path=NULL, format=c("pdf", "dvi"), clean=FALSE, quiet=TRUE, texinputs=NULL, ..., outPath=".", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- if (is.null(path)) filename else file.path(path, filename);
  if (!isUrl(pathname)) {
    pathname <- Arguments$getReadablePathname(pathname);
  }

  # Arguments 'outPath':
  outPath <- Arguments$getWritablePath(outPath);
  if (is.null(outPath)) outPath <- ".";

  # Arguments 'texinputs':
  texinputs <- Arguments$getCharacters(texinputs);

  # Arguments 'format':
  format <- match.arg(format);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Compiling LaTeX document");

  # Download URL?
  if (isUrl(pathname)) {
    verbose && enter(verbose, "Downloading URL");
    url <- pathname;
    verbose && cat(verbose, "URL: ", url);
    pathname <- downloadFile(url, verbose=less(verbose,50));
    verbose && cat(verbose, "Local file: ", pathname);
    verbose && exit(verbose);
  }

  ## Assert supported filename extension
  ext <- file_ext(pathname)
  if (!ext %in% c("tex", "ltx")) {
    throw("Unknown LaTeX filename extension (should lower case *.tex or *.ltx): ", pathname)
  }

  # Shorten, e.g. ../foo/../foo/ to ../foo
  pathname <- normalizePath(pathname);
  pathname <- getAbsolutePath(pathname);
  verbose && cat(verbose, "LaTeX pathname (absolute): ", pathname);
  verbose && printf(verbose, "Input file size: %g bytes\n", file.info(pathname)$size);
  verbose && cat(verbose, "Output format: ", format);
  verbose && cat(verbose, "Output and working directory: ", getAbsolutePath(outPath));
  pattern <- "(.*)[.]([^.]+)$";
  replace <- sprintf("\\1.%s", format);
  filenameOut <- gsub(pattern, replace, basename(pathname));
  pathnameOut <- filePath(outPath, filenameOut);
  verbose && cat(verbose, "Output pathname (", toupper(format), "): ", getAbsolutePath(pathnameOut));

  opwd <- ".";
  on.exit(setwd(opwd), add=TRUE);
  if (!is.null(outPath)) {
    opwd <- setwd(outPath);
  }

  verbose && enter(verbose, "Calling tools::texi2dvi()");
  pdf <- (format == "pdf");
  pathnameR <- getRelativePath(pathname);
  # Sanity check
  pathnameRx <- Arguments$getReadablePathname(pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Append the directory of the TeX file to TEXINPUTS search path?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathR <- dirname(pathnameR);
  if (pathR != ".") {
    verbose && enter(verbose, "Appending directory of TeX file to 'texinputs'");
    if (!is.null(texinputs)) {
      texinputs <- unlist(strsplit(texinputs, split="[:;]", fixed=FALSE), use.names=FALSE);
    }
    verbose && cat(verbose, "'texinputs' before:");
    verbose && print(verbose, texinputs);
    # Shorten, e.g. ../foo/../foo/ to ../foo
    pathR <- normalizePath(pathR);

    # Append as relative or absolute path (using the shortest one)
    pathRR <- getRelativePath(pathR);
    pathRA <- getAbsolutePath(pathR);
    if (nchar(pathRA) < nchar(pathRR)) {
      texinputs <- c(pathRA, texinputs);
    } else {
      texinputs <- c(pathRR, texinputs);
    }
##      texinputs <- c(pathRR, texinputs);
##      texinputs <- c(pathRA, texinputs);

    # Append as temporary link, iff possible
    verbose && enter(verbose, "Appending temporary link, iff possible");
    link <- basename(tempdir());
    verbose && cat(verbose, "Link: ", link);
    if (!file.exists(link)) {
      verbose && cat(verbose, "Trying to create link to target: ", pathR);
      tryCatch({
        linkT <- createLink(target=pathR, link=link);
        verbose && cat(verbose, "Created link: ", linkT);
      }, error = function(ex) {
        verbose && print(verbose, ex)
      })
      if (file.exists(link)) {
        linkA <- getAbsolutePath(link);
        on.exit(removeDirectory(linkA, mustExist=FALSE), add=TRUE);
        verbose && cat(verbose, "Link created: ", link);
        texinputs <- c(link, texinputs);
      }
    }
    verbose && exit(verbose);

    # Keep unique
    texinputs <- unique(texinputs);

    verbose && exit(verbose);
  }

  verbose && cat(verbose, "texinputs:");
  verbose && print(verbose, texinputs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Temporarily cleanup TEXINPUTS, BIBINPUTS, BSTINPUTS, TEXINDY
  # by removing empty, duplicated and non-existing paths.  This
  # lowers the risk for compilation failure due to too long paths.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Cleaning up LaTeX environment variable");
  cleanupPath <- function(path, sep=.Platform$path.sep) {
    appendSep <- (regexpr(sprintf("%s$", sep), path) != -1L);

    path <- unlist(strsplit(path, split=sep, fixed=TRUE));

    # Drop duplicates
    path <- unique(path);
    # Drop empty paths
    path <- path[nchar(path) > 0L];
    # Drop non-existing paths (accounting for foo// specified paths)
    pathX <- gsub("[/\\]*$", "", path);
    isDir <- sapply(pathX, FUN=file.exists);
    path <- path[isDir];
    # Re-append separator to the end?
    if (appendSep) path <- c(path, "");

    paste(path, collapse=sep);
  } # cleanupPath()

  vars <- c("TEXINPUTS", "BIBINPUTS", "BSTINPUTS", "TEXINDY");
  verbose && cat(verbose, "Original:");
  verbose && printf(verbose, " %s: %s\n", vars, Sys.getenv(vars));
  envs <- Sys.getenv(vars, NA);
  envs <- envs[!is.na(envs)];

  # Cleanup paths
  if (length(envs) > 0L) {
    # Undo any changes to system environments
    on.exit(do.call(Sys.setenv, as.list(envs)), add=TRUE);
    envs2 <- sapply(envs, FUN=cleanupPath);
    if (length(envs2) > 0L) {
      do.call(Sys.setenv, as.list(envs2))
      vars <- names(envs2)
      verbose && cat(verbose, "Cleaned up:");
      verbose && printf(verbose, " %s: %s\n", vars, Sys.getenv(vars));
    }
  }

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fake output PDF/DVI in case LaTeX is not available?
  # This makes it possible to test the compilation of vignettes
  # in 'R CMD check' up to the point of compiling a LaTeX file
  # into a PDF.  This is useful on Travis CI where we then can
  # avoid having to install a huge LaTeX system for each check.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fallback <- Sys.getenv("R_RSP_COMPILELATEX_FALLBACK");
  if (nzchar(fallback)) {
    verbose && enter(verbose, "Fallback");
    # Disable fallback until done to avoid recursive calls when
    # calling isCapableOf(..., "latex").
    Sys.unsetenv("R_RSP_COMPILELATEX_FALLBACK");
    fallback0 <- fallback;
    on.exit(Sys.setenv("R_RSP_COMPILELATEX_FALLBACK"=fallback0), add=TRUE);

    verbose && cat(verbose, "R_RSP_COMPILELATEX_FALLBACK=", fallback);
    forceFB <- (regexpr("-force", fallback) != -1L);
    if (forceFB) {
      fallback <- gsub("-force", "", fallback);
      verbose && cat(verbose, "Forced fallback");
      verbose && cat(verbose, "R_RSP_COMPILELATEX_FALLBACK=", fallback);
    }

    if (forceFB || !isCapableOf(R.rsp, "latex")) {
      if (fallback == "copy") {
        texi2dvi <- function(pathnameR, pdf=TRUE, ...) {
          verbose && enter(verbose, "Faking texi2dvi() by copying source file.");
          pathnameD <- file_path_sans_ext(basename(pathnameR));
          pathnameD <- paste(pathnameD, if (pdf) "pdf" else "dvi", sep=".");
          copyFile(pathnameR, pathnameD, overwrite=TRUE);
          verbose && printf(verbose, "Copied: %s -> %s\n", pathnameR, pathnameD);
          verbose && exit(verbose);
        } # texi2dvi()
        verbose && cat(verbose, "Faking texi2dvi() by copying source file.");
      } else {
        throw("Unknown value on _RSP_COMPILELATEX_FALLBACK_: ", fallback);
      }
    }
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Compile
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  texi2dvi(pathnameR, pdf=pdf, clean=clean, quiet=quiet, texinputs=texinputs);


  verbose && exit(verbose);

  setwd(opwd); opwd <- ".";
  verbose && printf(verbose, "Output file size: %g bytes\n", file.info(pathnameOut)$size);

  verbose && exit(verbose);

  pathnameOut;
}) # compileLaTeX()


############################################################################
# HISTORY:
# 2014-05-17
# o If environment variable _RSP_COMPILELATEX_FALLBACK_=copy and LaTeX
#   is not available, then compileLaTeX() will emulate the LaTeX compiler
#   by copying the input TeX file "as is" to the output PDF/DVI file,
#   e.g. foo.tex -> foo.pdf.
# 2014-04-06
# o ROBUSTNESS: Now compileLaTeX() cleans up and shortens LaTeX
#   environment variable iff possible before compiling the document.
# 2014-03-24
# o ROBUSTNESS: Now compileLaTeX() tries to shorten any paths as far
#   as possible, e.g. ../foo/../foo/ to ../foo/ to workaround possible
#   length limits of the TeX compiler.  It now also adds a symbolic link
#   to TEXINPUTS that refers to the directory of the LaTeX file.
# 2014-01-13
# o ROBUSTNESS: Now compileLaTeX() adds the directory of the LaTeX file
#   to TEXINPUTS also by its relative path (in addition to its absolute
#   path).  This provides a workaround for systems that does not handle
#   TEXINPUTS paths that are too long.  How to know what "too long" is is
#   not clear, but for the record a path with 138 characters is too long.
# 2013-07-16
# o Now compileLaTeX() adds the directory of the LaTeX file to the
#   TEXINPUTS search path, iff it's different than the working
#   directory.
# 2013-02-18
# o Added argument 'fake' to compileLaTeX().
# 2012-12-06
# o Added argument 'outPath' to compileLaTeX(), which is also the
#   working directory.
# o BUG FIX: compileLaTeX() would return an incorrect pathname unless
#   the given LaTeX file was in the current working directory.
# 2011-04-19
# o Added arguments 'clean' and 'quiet' to compileLaTeX().
# 2011-04-12
# o Added compileLaTeX().
# o Created.
############################################################################
