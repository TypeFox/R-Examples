###########################################################################/**
# @RdocFunction parseVignette
#
# @title "Parses an Rnw file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathname}{The Rnw file to be parsed.}
#   \item{commentPrefix}{A regular expression specifying the prefix
#     pattern of vignette comments.}
#   \item{final}{If @TRUE, the output PDF or HTML file is also located.}
#   \item{source}{If @TRUE, the output R source code file is also located.}
#   \item{maxLines}{The maximum number of lines to scan.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list or NULL if a non-vignette.
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
parseVignette <- function(pathname, commentPrefix="^[ \t]*%[ \t]*", final=FALSE, source=FALSE, maxLines=-1L, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  findOutput <- function(pathname, pattern) {
     path <- dirname(pathname);
     filename <- basename(pathname);
     ext <- gsub(".*[.]([^.]*)$", "\\1", filename);

     # All available output files
     filenames <- list.files(path=path, pattern=pattern);
     fullnames <- gsub("[.][^.]*$", "", filenames);
     patterns <- sprintf("^%s.*[.]%s$", fullnames, ext);
     keep <- (unlist(lapply(patterns, FUN=regexpr, filename), use.names=FALSE) != -1L);
     filenames <- filenames[keep];
     if (length(filenames) == 0L) return(NULL);

     # Order by decreasing filename lengths
     o <- order(nchar(filenames), decreasing=TRUE);
     filenames <- filenames[o];
     file.path(path, filenames);
  }


  if (!file.exists(pathname)) {
    stop("Cannot build vignette. File not found: ", pathname);
  }

  bfr <- readLines(pathname, warn=FALSE, n=maxLines);

  # Parse for "\Vignette" options
  pattern <- sprintf("%s\\\\Vignette(.*)\\{(.*)\\}", commentPrefix);
  rows <- which(regexpr(pattern, bfr) != -1L);
  bfr <- bfr[rows];

  # Nothing found?
  if (length(bfr) == 0L) {
    return(NULL);
  }

  # If the first entry is not among the first 20 rows, assume the ones
  # founds are part of the text document such entries rather than entries
  # used for the vignette itself.
  if (rows[1L] > 20L) {
    return(NULL);
  }

  opts <- grep(pattern, bfr, value=TRUE);
  keys <- gsub(pattern, "\\1", opts);
  values <- gsub(pattern, "\\2", opts);
  names(values) <- keys;
  opts <- as.list(values);

  # Drop duplicated entries, assuming the first ones are the intended
  # ones.  The extra ones may happen when a vignette documents how to
  # use %\\VignetteNnn{} markup.
  keep <- !duplicated(names(opts));
  opts <- opts[keep];

  # No %\VignetteIndexEntry{}?
  if (!is.element("IndexEntry", names(values))) {
    return(NULL);
  }

  vign <- c(list(pathname=pathname), opts);


  # Look for a generated PDF or HTML file?
  if (final) {
     output <- findOutput(pathname, pattern="[.](pdf|PDF|html|HTML)$");
     if (length(output) == 0L) {
       stop("Failed to located PDF or HTML output file for vignette: ", pathname);
     } else if (length(output) > 1L) {
       stop("Located more than one PDF or HTML output file for vignette: ", pathname);
     }
     vign$final <- output;
  }

  # Look for a generated R source code file?
  if (source) {
     output <- findOutput(pathname, pattern="[.][rRsS]$");
     if (length(output) > 1L) {
       output <- output[1L];
     }
     vign$source <- output;
  }

  # Assert unique entries
  names <- names(vign);
  dups <- names[duplicated(names)];
  if (length(dups) > 0L) {
    throw("Duplicated entries detected: ", paste(dups, collapse=", "));
  }

  vign;
} # parseVignette()



###########################################################################/**
# @RdocFunction parseVignettes
#
# @title "Locates and parses all vignettes"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The directory where to search for vignettes.}
#   \item{pattern}{Filename pattern to locate vignettes.}
#   \item{...}{Additional arguments passed to @see "parseVignette".}
#   \item{drop}{A @vector of filename patterns of vignette sources
#    to be ignored.}
# }
#
# \value{
#   Returns a @list where each element corresponds to an
#   identified vignette source file.  A file is considered to be
#   a vignette source file if it has \code{\\Vignette.*\{\}} markups
#   in the top 50 lines.
#   Each such vignette element consists of a named @list with
#   the parse \code{\\Vignette.*\{\}} information.
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
parseVignettes <- function(path=".", pattern="[.][^.~]*$", ..., drop="^dummy.tex$") {
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE);

  # Ignore certain files, e.g. "^dummy.Rnw$"?
  if (length(drop) > 0L) {
    filenames <- basename(pathnames);
    excl <- rep(FALSE, times=length(filenames));
    for (pattern in drop) {
      excl <- excl | (regexpr(pattern, filenames) != -1L);
    }
    pathnames <- pathnames[!excl];
  }

  vigns <- list();
  for (kk in seq_along(pathnames)) {
    pathname <- pathnames[kk];
    vign <- parseVignette(pathname, ...);
    if (length(vign) == 0L)
       next;
    vigns[[pathname]] <- vign;
  }

  vigns;
} # parseVignettes()



###########################################################################/**
# @RdocFunction buildNonSweaveVignette
#
# @title "Builds a non-Sweave Rnw vignette"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathname}{The vignette file to be built.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) a named @list.
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
buildNonSweaveVignette <- function(vign, envir=new.env(), ...) {
  # Local functions
  SweaveStangle <- function(file, ...) {
    pathnameR <- utils::Sweave(file, ...);
    utils::Stangle(file, ...);
    pathnameR;
  }

  # A filename?
  if (is.character(vign)) {
    pathname <- vign;
    vign <- parseVignette(pathname, ...);
  }

  pathname <- vign$pathname;

  # Load required packages
  if (!is.null(vign$Depends)) {
    pkgNames <- vign$Depends;
    pkgNames <- unlist(strsplit(pkgNames, split=",", fixed=TRUE), use.names=FALSE);
    pkgNames <- gsub("(^[ \t]*|[ \t]*$)", "", pkgNames);
    for (pkgName in pkgNames) {
      library(pkgName, character.only=TRUE);
    }
  }

  # Build vignette according to \VignetteBuild{} command
  if (!is.null(cmd <- vign$Engine) && nchar(cmd) > 0L) {
    # Retrieve the "engine" according to \VignetteEngine{} expression
    res <- get(cmd, envir=envir, mode="function");
  } else if (!is.null(cmd <- vign$Build) && nchar(cmd) > 0L) {
    # Parse \VignetteBuild{} expression
    tryCatch({
      expr <- parse(text=cmd);
    }, error = function(ex) {
      stop(sprintf("Syntax error in \\VignetteBuild{%s}: %s", cmd, ex$message));
    });

    # Evaluate \VignetteBuild{} expression
    res <- eval(expr);
  } else {
     # If not specified, assume Sweave
     res <- SweaveStangle;
  }

  # Was a function specified?
  if (is.function(res)) {
    fcn <- res;
    res <- fcn(pathname);
  }

  invisible(res);
} # buildNonSweaveVignette()



###########################################################################/**
# @RdocFunction buildNonSweaveVignettes
#
# @title "Builds all non-Sweave Rnw vignette"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The directory where to search for non-Sweave vignettes.}
#   \item{pattern}{Filename pattern to locate non-Sweave vignettes.}
#   \item{...}{Additional arguments passed to @see "buildNonSweaveVignette".}
# }
#
# \value{
#   Returns (invisibly) a named @list with elements of what
#   the vignette builder returns.
# }
#
# @author
#
# \seealso{
#   To build vignette, see @see "buildNonSweaveVignette".
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
buildNonSweaveVignettes <- function(...) {
  vigns <- parseVignettes(...);
  if (length(vigns) > 0L) {
     envir <- new.env();
     path <- system.file("doc", "templates", package="R.rsp");
     path <- c(path, dirname(vigns[[1L]]$pathname));
     pathnames <- file.path(path, "enginesMap.R");
     pathnames <- pathnames[file_test("-f", pathnames)];
     for (pathname in pathnames) {
       expr <- parse(pathname);
       eval(expr, envir=envir);
     }
  }
  for (kk in seq_along(vigns)) {
    vign <- vigns[[kk]];
    vign$result <- buildNonSweaveVignette(vign, envir=envir, ...);
    vigns[[kk]] <- vign;
  }
  invisible(vigns);
} # buildNonSweaveVignettes()



###########################################################################/**
# @RdocFunction buildNonSweaveTexToPdf
#
# @title "Compiles all TeX files into PDFs"
#
# \description{
#  @get "title", unless already done.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The directory where to search for TeX files.}
#   \item{pattern}{Filename pattern to locate TeX files.}
#   \item{...}{Additional arguments passed to @see "tools::texi2pdf".}
# }
#
# \value{
#   Returns (invisibly) a named @list of results.
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
buildNonSweaveTexToPdf <- function(path=".", pattern="[.](tex|ltx)$", ...) {
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE);

  # Ignore dummy.tex (which is created by R from dummy.Rnw just before make)
  keep <- !is.element(basename(pathnames), c("dummy.tex"));
  pathnames <- pathnames[keep];

  res <- list();
  for (pathname in pathnames) {
    pathnamePDF <- sprintf("%s.pdf", gsub(pattern, "", pathname));
    if (!isFile(pathnamePDF)) {
       res[[pathname]] <- texi2pdf(file=pathname, ...);
    }
  }
  invisible(res);
} # buildNonSweaveTexToPdf()


###########################################################################/**
# @RdocFunction buildPkgIndexHtml
#
# @title "Builds a package index HTML file"
#
# \description{
#  @get "title", iff missing.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the absolute pathame to the built index.html file.
#   If an index.html file already exists, nothing is done and @NULL
#   is returned.
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
buildPkgIndexHtml <- function(...) {
  # Nothing to do?
  if (file.exists("index.html")) {
    return(NULL);
  }

  library("R.rsp");

  filename <- "index.html.rsp";
  if (!file.exists(filename)) {
    # If not custom index.html.rsp exists, use the one of the R.rsp package
    path <- system.file("doc/templates", package="R.rsp");
    pathname <- file.path(path, filename);
    file.copy(pathname, to=".");
    on.exit({
      file.remove(filename);
    });
  }

  # Sanity check
  stopifnot(file.exists(filename));

  # Build index.html
  rfile(filename);
} # buildPkgIndexHtml()


############################################################################
# HISTORY:
# 2014-04-30
# o ROBUSTNESS: Now parseVignette() drops duplicated %\VignetteNnn{}
#   entries, extra entries that may appear because the vignette actually
#   talks about such markup, like some of the vignette in this package.
# 2013-10-13
# o BUG FIX: parseVignette() ignores files that do not contain a
#   '%\VignetteIndexEntry{}'.
# o Added argument 'dropDummy' to parseVignettes().
# 2013-03-28
# o Now buildNonSweaveTexToPdf() ignores 'dummy.tex'.
# 2013-03-07
# o Deprecated use of \VignetteBuild{} in favor of \VignetteEngine{}
#   together with an 'enginesMap.R' file.
# o Dropped use of \VignetteSource{}.
# o Added parseVignettes().
# o Now parseVignette() only scans the first 50 lines.
# 2013-02-14
# o Added Rdoc help for all functions.
# o Now buildNonSweaveVignette() also handles \VignetteBuild{R.rsp::rfile}
#   given that \VignetteSource{} is specified.
# 2011-11-23
# o Added buildPkgIndexHtml().
# o Added parseVignette().
# o Added buildNonSweaveVignettes().
# o Added buildNonSweaveVignette().
# o Created.
############################################################################
