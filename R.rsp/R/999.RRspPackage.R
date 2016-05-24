###########################################################################/**
# @RdocClass RRspPackage
#
# @title "The RRspPackage class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setConstructorS3("RRspPackage", function(...) {
  extend(Package(...), "RRspPackage");
})



###########################################################################/**
# @RdocMethod capabilitiesOf
# @aliasmethod isCapableOf
#
# @title "Checks which tools are supported"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{Optional @character @vector of which tools to check.}
#  \item{force}{If @TRUE, cached results are ignored, otherwise not.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @logical named @character @vector.
# }
#
# \examples{
#   # Display which tools are supported by the package
#   print(capabilitiesOf(R.rsp))
#
#   # Check whether AsciiDoc is supported
#   print(isCapableOf(R.rsp, "asciidoc"))
#
#   # Check whether pandoc v1.12 or newer is supported
#   print(isCapableOf(R.rsp, "pandoc (>= 1.12)"))
# }
#
# @author "HB"
#
#*/###########################################################################
setMethodS3("capabilitiesOf", "RRspPackage", function(static, what=NULL, force=FALSE, ...) {
  res <- static$.capabilities;
  if (force || is.null(res)) {
    res <- list();

    # Check software
    res$asciidoc <- !is.null(findAsciiDoc(mustExist=FALSE));
    res$knitr <- !is.null(isPackageInstalled("knitr"));
    res$markdown <- !is.null(isPackageInstalled("markdown"));
    res$pandoc <- !is.null(findPandoc(mustExist=FALSE));
    res$sweave <- !is.null(isPackageInstalled("utils"));

    # Check LaTeX
    path <- system.file("rsp_LoremIpsum", package="R.rsp");
    pathname <- file.path(path, "LoremIpsum.tex");
    res$latex <- tryCatch({
      pathnameR <- compileLaTeX(pathname, outPath=tempdir());
      isFile(pathnameR);
    }, error = function(ex) FALSE);

    # Order lexicographically
    o <- order(names(res));
    res <- res[o];

    # Coerce into a named character vector
    res <- unlist(res, use.names=TRUE);

    # Record
    static$.capabilities <- res;
  }

  if (!is.null(what)) {
    res <- res[what];
  }

  res;
}, static=TRUE)


setMethodS3("isCapableOf", "RRspPackage", function(static, what, ...) {
  # Argument 'what':
  what <- Arguments$getCharacter(what);
  pattern <- "^([^ ]+)[ ]*(|[(](<|<=|==|>=|>)[ ]*([^)]+)[)])$";
  if (regexpr(pattern, what) == -1L) {
    throw("Unknown syntax of argument 'what': ", what);
  }

  name <- gsub(pattern, "\\1", what);
  op <- gsub(pattern, "\\3", what);
  ver <- gsub(pattern, "\\4", what);
  if (nzchar(op)) {
    op <- get(op, mode="function", envir=baseenv());
  } else if (nzchar(ver)) {
    throw("Missing version operator in argument 'what': ", what);
  }

  res <- capabilitiesOf(static, what=name, ...);

  # Nothing more to do?
  if (!is.element(name, names(res))) {
    return(FALSE);
  }

  # Nothing more to do?
  if (!nzchar(ver)) {
    return(res);
  }

  # Get available version
  if (name == "asciidoc") {
    v <- attr(findAsciiDoc(mustExist=FALSE), "version");
  } else if (name == "pandoc") {
    v <- attr(findPandoc(mustExist=FALSE), "version");
  } else if (is.element(name, c("knitr", "markdown"))) {
    v <- packageVersion(name);
  } else if (name == "sweave") {
    v <- packageVersion("utils");
  } else {
    v <- NA;
  }

  # Compare to requested version
  res <- isTRUE(op(v, ver));

  res;
})


############################################################################
# HISTORY:
# 2013-09-28
# o No isCapableOf() also supports version specifications.
#   Ideally capabilitiesOf() and isCapableOf() should be moved to
#   the Package class of the 'R.oo' package.
# 2013-07-19
# o Created from AromaSeq.R in aroma.seq.
############################################################################
