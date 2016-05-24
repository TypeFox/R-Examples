##########################################################################/**
# @RdocClass Rdoc
#
# @title "Class for converting Rdoc comments to Rd files"
#
# \description{
#  @classhierarchy
#
#  @get "title".
# }
#
# @synopsis
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# \examples{\dontrun{@include "../incl/Rdoc.Rex"}}
#
# \references{
#   R developers,
#   \emph{Guidelines for Rd files},
#   \url{http://developer.r-project.org/Rds.html},
#   2003
# }
#
# @keyword documentation
#*/###########################################################################
setConstructorS3("Rdoc", function() {
  extend(Object(), "Rdoc",
    .START      = paste(sep="", "/", "**"),  # To hide it from itself!!!
    .STOP       = paste(sep="", "*", "/"),   #   -  "  -
    .nameFormat = "method.class",
    .manPath    = "../man/",
    package     = NULL,
    source      = NULL
  )
})



###########################################################################/**
# @RdocMethod getNameFormat
#
# @title "Gets the current name format"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "setNameFormat"
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getNameFormat", "Rdoc", function(static, ...) {
  Rdoc$.nameFormat;
}, static=TRUE)




###########################################################################/**
# @RdocMethod setNameFormat
#
# @title "Sets the current name format"
#
# \description{
#   @get "title".
#   Throws a @see "RccViolationException" if an unknown format is requested.
# }
#
# @synopsis
#
# \arguments{
#   \item{nameFormat}{
#    If \code{"method.class"}, help files for methods belonging to classes are
#    named <method>.<class>.
#    If \code{"class.method"}, help files for methods belonging to classes are
#    named <class>.<method>.
#    These are currently the only name formats supported.
#   }
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
#   @seemethod "getNameFormat"
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("setNameFormat", "Rdoc", function(static, nameFormat, ...) {
  if (nameFormat == "class.method") {
  } else if (nameFormat == "method.class") {
  } else {
    throw(RdocException("Unknown name format: ", nameFormat));
  }
  Rdoc$.nameFormat <- nameFormat;
}, static=TRUE)




###########################################################################/**
# @RdocMethod getKeywords
#
# @title "Gets the keywords defined in R with descriptions"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "setManPath"
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getKeywords", "Rdoc", function(this, fullInfo=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate the KEYWORDS.db file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- Sys.getenv("R_DOC_DIR");
  if (is.null(path) || nchar(path) == 0) {
    # Backward compatibility
    path <- file.path(Sys.getenv("R_HOME"), "doc");
    tryCatch({
      path <- R.home("doc");
    }, error = function(ex) {});

    if (!file.exists(path)) {
      throw("Cannot determine the R doc directory. R_DOC_DIR was not set and R_HOME/doc/ does not exist: ", path);
    }
  }

  pathname <- file.path(path, "KEYWORDS.db");
  if (!file.exists(pathname)) {
    throw("The KEYWORDS.db file was not found: ", pathname);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read keywords
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keywords <- readLines(pathname, warn=FALSE);
  keywords <- strsplit(keywords, ":");
  names <- lapply(keywords, FUN=function(x) x[1]);
  names <- unlist(names);
  names <- gsub("[[:space:]]+$", "", gsub("^[[:space:]]+", "", names));
  desc  <- lapply(keywords, FUN=function(x) x[2]);
  desc <- unlist(desc);
  desc <- gsub("[[:space:]]+$", "", gsub("^[[:space:]]+", "", desc));
  keywords <- strsplit(names, "\\|");
  len <- unlist(lapply(keywords, FUN=length));
  keywords <- unlist(lapply(keywords, FUN=function(x) x[length(x)]));
  keywords <- keywords[len > 1];
  desc <- desc[len > 1];
  names(keywords) <- desc;

  keywords;
}, static=TRUE);


###########################################################################/**
# @RdocMethod isKeyword
#
# @title "Checks if a word is a Rd keyword"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @logical.
# }
#
# @author
#
# \seealso{
#   @seemethod "getKeywords"
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("isKeyword", "Rdoc", function(this, word, ...) {
  is.element(word, Rdoc$getKeywords());
}, static=TRUE)



###########################################################################/**
# @RdocMethod getManPath
#
# @title "Gets the path to the directory where the Rd files will be saved"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "setManPath"
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getManPath", "Rdoc", function(this, ...) {
  this$.manPath;
}, static=TRUE);



###########################################################################/**
# @RdocMethod setManPath
#
# @title "Sets the path to the directory where the Rd files should be saved"
#
# \description{
#   @get "title". By default the path is \code{../man/} assuming that the
#   current directory is \code{../R/}, which is where source files commonly
#   are placed.
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
#   @seemethod "getManPath"
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("setManPath", "Rdoc", function(this, path="../man/", ...) {
  info <- file.info(path);
  if (is.na(info$isdir))
    path <- gsub("/$", "", path);
  info <- file.info(path);
  if (is.na(info$isdir))
    throw("Path does not exists: ", path);
  if (info$isdir != TRUE)
    throw("Specified path is not a directory: ", path);
  this$.manPath <- as.character(path);
}, static=TRUE);



###########################################################################/**
# @RdocMethod createManPath
#
# @title "Creates the directory where the Rd files should be saved"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns @TRUE if the directory was creates, @FALSE if it already exists
#  and throws an @Exception if failed.
# }
#
# @author
#
# \seealso{
#   @seemethod "getManPath"
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("createManPath", "Rdoc", function(this, ...) {
  # Check if the path already exists, otherwise create it.
  path <- getManPath(this);
  # file.exists() and file.info() is sensitive to trailing '/'.
  path <- gsub("/$", "", path);
  isdir <- as.logical(file.info(path)["isdir"]);
  if (file.exists(path) && isdir)
    return(FALSE);

  # Path 'destPath' does not exist at all. It might be that there is
  # a file with the same name, but in any case, the OS should
  # take care of conflict if it thinks it is a conflict.
  if (!dir.create(path))
    throw(Exception("Could not create destination directory: ", path));

  return(TRUE);
}, static=TRUE, protected=TRUE)





###########################################################################/**
# @RdocMethod createName
#
# @title "Creates a class-method name"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{class}{A class name (@character string).}
#   \item{method}{A method name (@character string).}
#   \item{escape}{If @TRUE, non-valid filename characters are escaped into
#     valid character strings.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "escapeRdFilename".
#   @seeclass
# }
#
# @keyword documentation
# @keyword internal
#*/###########################################################################
setMethodS3("createName", "Rdoc", function(static, class, method, escape=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get name format to be used (can be set globally)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nameFormat <- Rdoc$getNameFormat();
  if (nameFormat == "class.method") {
    name <- paste(class, ".", method, sep="");
  } else if (nameFormat == "method.class") {
    name <- paste(method, ".", class, sep="");
  } else {
    throw(RdocException("Unknown name format: ", nameFormat));
  }

  if (escape) {
    name <- Rdoc$escapeRdFilename(name);
  }

  name;
}, static=TRUE, private=TRUE);



###########################################################################/**
# @RdocMethod escapeRdFilename
#
# @title "Escape non-valid characters in a filename"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{A filename (@character string) to be escaped.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("escapeRdFilename", "Rdoc", function(static, filename, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Escape non-valid filenames
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- gsub("\\\\\\$", "DOLLAR", filename);
  filename <- gsub("[$]", "DOLLAR", filename);
  filename <- gsub("<-", "< -", filename);
  filename <- gsub("<", "LT", filename);
  filename <- gsub("[[]", "-LB-", filename);

  # From R v1.8.1 can't Rd filenames contain whitespace.
  filename <- gsub("[ \t]", "_", filename);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # A filename must start with a letter or a digit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  startOk <- (regexpr("^[abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]", filename) != -1);
  if (!startOk) {
    # Fix Rd filename (not really important actually).
    filename <- paste("000", filename, sep="");
  }

  filename;
}, protected=TRUE, static=TRUE) # escapeRdFilename()




###########################################################################/**
# @RdocMethod compile
#
# @title "Compile source code files containing Rdoc comments into Rd files"
#
# \description{
#   @get "title".
#
#  \emph{Note, the class and methods to be compiled have to be loaded into
#  \R by for instance \code{library()} or  \code{source()} before calling
#  this method.}
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The pathname or filename pattern of the Rdoc files to be
#     compiled.}
#   \item{destPath}{The path where the generated Rd files should be saved.}
#   \item{showDeprecated}{If @TRUE, Rd files are generated for deprecated
#     objects too, otherwise not.}
#   \item{addTimestamp}{If @TRUE, a date and time stamp is added to the
#     Rd header comments.  This timestamp might be confusing for version
#     control systems, which is why it can be turned off with @FALSE.}
#   \item{locale}{The locale to be set/used when compiling Rdoc comments.
#     This help assuring strings are sorted the same way across systems.}
#   \item{source}{If @TRUE, the Rdoc files will be \code{source()}:ed first.
#     This work of course only for Rdoc files that are R source files.}
#   \item{verbose}{If @TRUE, detailed compilation information is printed.}
#   \item{debug}{If @TRUE, extra debug information is printed.}
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
# @keyword documentation
#*/###########################################################################
setMethodS3("compile", "Rdoc", function(this, filename=".*[.]R$", destPath=getManPath(this), showDeprecated=FALSE, addTimestamp=FALSE, locale="C", verbose=FALSE, source=FALSE, check=TRUE, debug=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Global variables
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  authorWarn <- FALSE;
  pkgAuthors <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isCapitalized <- function(str) {
    first <- substring(str,1,1);
    (first == toupper(first))
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # a d d K e y w o r d ( )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  validateKeyword <- function(keyword) {
    knownKeywords <- Rdoc$getKeywords();
    if (!is.element(keyword, knownKeywords)) {
      alts <- agrep(keyword, knownKeywords);
      alts <- paste("'", knownKeywords[alts], "'", collapse=", ", sep="");
      if (nchar(alts) > 0)
        alts <- paste("Did you mean ", alts, "?", sep="");
      throw(RdocException("Unknown keyword: ", keyword, ". ",
                                                alts, source=sourcefile));
    }
  } # validateKeyword()


  rdocKeywords <- c();

  addKeyword <- function(keyword) {
    keyword <- as.character(keyword);

    # A remove keyword?
    if (regexpr("^-", keyword) != -1) {
      rdocKeywords <<- unique(c(rdocKeywords, keyword));
      keyword <- gsub("^-", "", keyword);
    } else {
      rdocKeywords <<- unique(c(rdocKeywords, keyword));
    }

    # Validate keyword
    validateKeyword(keyword);
  } # addKeyword()


  getRdKeywords <- function(...) {
    # Get all keywords
    if (length(rdocKeywords) == 0)
      return("");

    isRemove <- (regexpr("^-", rdocKeywords) != -1);
    keywords <- rdocKeywords[!isRemove];
    exclKeywords <- gsub("^-", "", rdocKeywords[isRemove]);
    keywords <- setdiff(keywords, exclKeywords);
    keywords <- unique(keywords);

    # Empty current list of keywords
    rdocKeywords <<- c();

    if (length(keywords) == 0)
      return(NULL);

    lines <- paste("\\keyword{", keywords, "}", sep="");
    lines <- paste(lines, collapse="\n");
    lines;
  } # getRdKeywords()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # e s c a p e N a m e ( )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # From the help:
  #  \name{name}
  #    name typically is the basename of the Rd file containing the
  #    documentation. It is the "name" of the Rd object represented
  #    by the file, has to be unique in a package, and must not
  #    contain LaTeX special characters (#, $, %, &, ~, _, ^, \, {, }).
  escapeName <- function(name) {
    name <- gsub("\\#", "POUND", name);
    name <- gsub("\\$", "DOLLAR", name);
    name <- gsub("\\%", "PERCENT", name);
    name <- gsub("\\&", "AND", name);
    name <- gsub("\\~", "TILDE", name);
    name <- gsub("\\_", "UNDERSCORE", name);
    name <- gsub("\\^", "POWER", name);
    name <- gsub("\\\\", "BACKSLASH", name);
    name <- gsub("\\{", "LCURLY", name);
    name <- gsub("\\}", "RCURLY", name);
    name <- gsub("<-", "< -", name);
    name;
  } # escapeName()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # e s c a p e A l i a s ( )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  escapeAlias <- function(alias) {
    # Don't escape aliases?!?  /HB 2004-03-03
    alias <- gsub("\\%", "\\\\%", alias);
#    alias <- gsub("\\$", "\\\\$", alias);
#    alias <- gsub("<-", "< -", alias);
    alias;
  } # escapeAlias()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # e x t r a c t R d o c s ( )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extractRdocs <- function(filename, verbose=FALSE, debug=FALSE) {
    if (!file.exists(filename))
      throw(RdocException("File not found: ", filename));

    # Read all lines from the source code file
    lines <- readLines(filename, warn=FALSE);

    if (length(lines) == 0)
      return(list());

    # Keep only the lines that are comments.
    lines <- lines[(regexpr("^ *#", lines) != -1)];

    if (length(lines) == 0)
      return(list());

    # Find all Rdoc begins and ends
    begins <- which(regexpr("/\\*\\*", lines) != -1);
    ends <- which(regexpr("\\*/", lines) != -1);

    if (length(begins) != length(ends))
      throw(RdocException("Number of Rdoc begins do not match number of Rdoc ends: ", filename));

    if (any(begins - ends > 0))
      throw(RdocException("Some of the Rdoc begins comes after the Rdoc ends: ", filename));


    rdocs <- list();
    for (k in seq_along(begins)) {
      idx <- begins[k]:ends[k];
      tmp <- lines[idx];

      # Remove everything before the begin tag including the tag
      tmp[1] <- gsub("^#.*/\\*\\*", "", tmp[1]);

      # Remove everything after the end tag including the tag
      last <- length(tmp);
      tmp[last] <- gsub("^#.*\\*/.*", "", tmp[last]);

      # Remove all leading single and double comment characters
      tmp <- gsub("^#{1,3}", "", tmp);

      # Find (minimum) indentation
      indents <- regexpr("[^ ]", tmp[nchar(tmp) > 0])
      indent <- min(indents);
      tmp <- substring(tmp, first=indent);

      # Remove all trailing whitespace
      tmp <- gsub("[ \t\v\r\n]$", "", tmp);

      # Remove all empty lines at the beginning
      while (nchar(tmp[1]) == 0)
        tmp <- tmp[-1];

      # Remove all empty lines at the end
      while (nchar(tmp[length(tmp)]) == 0)
        tmp <- tmp[-length(tmp)];

      attr(tmp, "sourcefile") <- filename;
      rdocs[[k]] <- tmp;
    }

    rdocs;
  } # extractRdocs()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # w r i t e R d ( )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  sourcefile <- NA;
  writeRd <- function(rds, path=getManPath(this), addTimestamp=TRUE, verbose=FALSE, debug=FALSE) {
    for (rd in rds) {
      name <- attr(rd, "name");
      if (!is.null(path)) {
        if (regexpr("/$", path) == -1 && regexpr("\\$", path) == -1)
          path <- paste(path, "/", sep="");
      }

      filename <- Rdoc$escapeRdFilename(name);
      filename <- paste(path, filename, ".Rd", sep="");
      if (verbose) {
        cat("Generating ", filename, "...", sep="");
      }
      sourcefile <<- sourcefile <- attr(rd, "sourcefile");

      hdr <- c("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
      hdr <- c(hdr, "% Do not modify this file since it was automatically generated from:");
      hdr <- c(hdr, "% ");
      hdr <- c(hdr, paste("%  ", sourcefile, sep=""));
      hdr <- c(hdr, "% ");
      if (addTimestamp) {
        hdr <- c(hdr, paste("% on ", date(), ".", sep=""));
        hdr <- c(hdr, "% ");
      }
      hdr <- c(hdr, "% by the Rdoc compiler part of the R.oo package.");
      hdr <- c(hdr, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
      bfr <- paste(c(hdr, "", rd, ""), collapse="\n");
      writeChar(bfr, eos=NULL, con=filename);
    }
  } # writeRd()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # c o m p i l e R d o c ()
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  compileRdoc <- function(rdocs, showDeprecated=FALSE, verbose=FALSE, debug=FALSE) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!is.list(rdocs))
      throw(RdocException("Internal error: Expected a list: ", class(rdocs)[1L]));


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Defines all simple tags that are shortcuts to different help documents.
    # Example: @TRUE -> \\code{\\link[base:logical]{TRUE}}
    # NA means that the text should just be place inside \code{}.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    shorttags <- list(
      # Package base:
      "0"           = NA,
      "1"           = NA,
      "..."         = NA,
      "."           = NA,         # Must come *after* "..."
      "array"       = "base",
      "attributes"  = "base",
      "attr"        = "base",
      "ANY"         = NA,
      "call"        = "base",
      "character"   = "base",
      "complex"     = "base",
      "connection"  = "base:connections",
      "data.frame"  = "base",
      "dimnames"    = "base",
      "dim"         = "base",
      "double"      = "base",
      "environment" = "base",
      "expression"  = "base",
      "factor"      = "base",
      "FALSE"       = "base:logical",
      "formula"     = "base",
      "function"    = "base",
      "Inf"         = "base:is.finite",
      "integer"     = "base",
      "length"      = "base",
      "list"        = "base",
      "logical"     = "base",
      "matrix"      = "base",
      "names"       = "base",
      "name"        = "base",
      "NA"          = "base",
      "NaN"         = "base:is.finite",
      "NULL"        = "base",
      "numeric"     = "base",
      "table"       = "base",
      "TRUE"        = "base:logical",
      "raw"         = "base",
      "ts"          = "base",
      "vector"      = "base",
      "warning"     = "base",
      # Package R.oo:
      "Object"      = "R.oo",
      "Exception"   = "R.oo",
      "throw"       = "R.oo"
    );

    names <- names(shorttags);
    match <- gsub("\\.", "\\\\.", names);
    attr(shorttags, "beginsWith") <- paste("^@", match, sep="");
    attr(shorttags, "contains") <- paste("[^abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9.]@", match, sep="");
    replace <- paste("\\link[", unlist(shorttags), "]{", names, "}", sep="");
    replace[is.na(shorttags)] <- names[is.na(shorttags)];
    replace <- paste("\\code{", replace, "}", sep="");
    attr(shorttags, "replace") <- replace;


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    getTagValue <- function(bfr) {
      # 1. Remove all leading whitespace
      bfr <- gsub("^[ \t]", "", bfr);

      # 2a. Is there a '{' + '}' pair?  (nesting brackets are not allowed)
      if ((beginPos <- regexpr("^\\{", bfr)) != -1L) {
        # Find the end '}'
        endPos <- regexpr("\\}", bfr);
        if (endPos == -1L)
          throw(RdocException("Closing } is missing: ", substring(bfr, first=1L, last=20L), source=sourcefile));
        value <- substring(bfr, first=beginPos+1L, last=endPos-1L);
        bfr <- substring(bfr, first=endPos+1L);
      }
      # 2b. ...or a '"' + '"' pair?  (*internal* escaped " are ignored)
      else if ((beginPos <- regexpr("^\"", bfr)) != -1L) {
        endPos <- regexpr("[^\\]\"", bfr);
        if (endPos == -1L)
          throw(RdocException("Closing \" is missing: ", substring(bfr, first=1L, last=20L), source=sourcefile));
        value <- substring(bfr, first=beginPos+1L, last=endPos);
        bfr <- substring(bfr, first=endPos+2L);
      }
      # 2c. ...otherwise the value is the first word found
      #     (on the same line!)
      else {
        beginPos <- 1L;
        endPos <- regexpr("([ \t\n\r]|$)", bfr);
        value <- substring(bfr, first=1L, last=endPos-1L);
        # Ad hoc. /HB 2013-03-25
        if (value != "}") {
          bfr <- substring(bfr, first=endPos);
        } else {
          value <- "";
        }
      }

      attr(bfr, "value") <- value;
      bfr;
    } # getTagValue()


    isObjectDeprecated <- function(name, ...) {
      obj <- getObject(this, name=name, ...);
      mods <- attr(obj, "modifiers");
      is.element("deprecated", mods);
    } # isObjectDeprecated()


    # Read and parse authors from DESCRIPTION's 'Authors@R' or 'Author'.
    getPackageAuthors <- function() {
      if (!is.null(pkgAuthors)) {
         return(pkgAuthors);
      }
      pkg <- Package(Rdoc$package);
      authors <- getAuthor(pkg, as="person");
      authorsN <- format(authors, include=c("given", "family"));

      maintainers <- getMaintainer(pkg, as="person");
      maintainersN <- format(maintainers, include=c("given", "family"));

      # Append maintainers, if not already listed as authors
      keep <- !is.element(maintainersN, authorsN);
      maintainers <- maintainers[keep];
      if (length(maintainers) > 0L) {
##        authors <- c(authors, maintainers);
      }

      authors;
    } # getPackageAuthors()


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocClass <- function(bfr) {
      bfr <- getTagValue(bfr);
      class <<- attr(bfr, "value");

      typeOfClass <- typeOfClass(class);
      if (is.na(typeOfClass)) {
        throw(RdocException("Class is either not defined or loaded, or not an S4/setClass() or S3/setConstructorS3() class: ", class));
      }

      if (typeOfClass == "S4") {
        clazz <<- getClass(class);
      } else if (typeOfClass == "S3-Object") {
        clazz <<- Class$forName(class);
      }

      line <- paste("\\name{", escapeName(class), "}\n", sep="");
      line <- paste(line, "\\docType{class}\n", sep="");
      line <- paste(line, "\\alias{", class, "}\n", sep="");
      addKeyword("classes");

      if (typeOfClass == "S4") {
        line <- paste(line, "\\alias{", class, "-class}\n", sep="");
        usage <- Rdoc$getClassS4Usage(clazz);
      } else if (typeOfClass == "S3-Object") {
        usage <<- Rdoc$getUsage(method=class);
        isDeprecated <<- isDeprecated(clazz);
      }

      name <<- class;
      objectName <<- class;

      rd <<- paste(rd, line, sep="");
      bfr;
    } # tagRdocClass()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagClass <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      class <<- value;
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocMethod <- function(bfr) {
      bfr <- getTagValue(bfr);
      method <- attr(bfr, "value");
      objectName <<- paste(method, class, sep=".");
      isDeprecated <<- isObjectDeprecated(objectName);

      # Find method
      fcn <- NULL;
      tryCatch({
        fcn <- Rdoc$getObject(objectName, mode="function");
      }, error = function(ex) {
        cat("Failed...\n");
        print(ex);
        cat("Failed...done\n");
      })

      if (!is.function(fcn)) {
        throw(RdocException("Could not get method. Function was not found: ", objectName, "()", source=Rdoc$source));
      }

      methodName <- createName.Rdoc(NULL, class, method, escape=FALSE);

      isStatic <- is.element("static", attr(fcn, "modifiers"));
      if (isStatic) {
        staticName <- paste(class, method, sep="$");
        name <- staticName;
        alias <- c(staticName, escapeAlias(methodName));
      } else {
        name <- escapeName(methodName);
        alias <- escapeAlias(methodName);
      }

      # Treat internal and non-internal methods differently
      if (isCapitalized(class)) {
        alias <- c(alias, paste(class, method, sep="."));
        alias <- c(alias, paste(method, ",", class, "-method", sep=""));
      }
      alias <- c(alias, paste(method, class, sep="."));

      # Multiple aliases(?)
      alias <- unique(alias);
      alias <- paste("\\alias{", alias, "}", sep="");

      line <- paste("\\name{", name, "}", sep="");
      line <- c(line, alias);
      line <- paste(line, collapse="\n");

      # Treat internal and non-internal methods differently
      if (isCapitalized(class)) {
        addKeyword("internal");
      }

      addKeyword("methods");

      name <<- methodName; # Filename
      usage <<- Rdoc$getUsage(method=method, class=class);
      rd <<- paste(rd, line, sep="");

      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocDefault <- function(bfr) {
      bfr <- getTagValue(bfr);
      default <- attr(bfr, "value");
      objectName <<- default;
      isDeprecated <<- isObjectDeprecated(objectName);
      name <- default;
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "}\n", sep="");
      line <- paste(line, "\\alias{", name, ".default}\n", sep="");
      line <- paste(line, "\\alias{", name, "}", sep="");
      usage <<- Rdoc$getUsage(method=default, class="default");
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocGeneric <- function(bfr) {
      bfr <- getTagValue(bfr);
      generic <- attr(bfr, "value");
      objectName <<- generic;
      isDeprecated <<- isObjectDeprecated(objectName);
      name <- generic;
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "}\n", sep="");
      line <- paste(line, "\\alias{", name, "}\n", sep="");
      addKeyword("methods");
      usage <<- Rdoc$getUsage(method=generic);
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocFunction <- function(bfr) {
      bfr <- getTagValue(bfr);
      fcn <- attr(bfr, "value");
      objectName <<- fcn;
      isDeprecated <<- isObjectDeprecated(objectName);
      name <- fcn;
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "}\n", sep="");
      line <- paste(line, "\\alias{", name, "}\n", sep="");
      usage <<- Rdoc$getUsage(method=fcn);
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocObject <- function(bfr) {
      bfr <- getTagValue(bfr);
      object <- attr(bfr, "value");
      objectName <<- object;
      isDeprecated <<- isObjectDeprecated(objectName);
      name <- object;
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "}\n", sep="");
      line <- paste(line, "\\alias{", name, "}", sep="");
      usage <<- NULL;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocData <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      name <- value;
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "}\n", sep="");
      line <- paste(line, "\\alias{", name, "}\n", sep="");
      line <- paste(line, "\\docType{data}\n", sep="");
      addKeyword("datasets");
      objectName <<- value;
      usage <<- NULL;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocDocumentation <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      name <- value;
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "}\n", sep="");
      line <- paste(line, "\\alias{", name, "}\n", sep="");
      addKeyword("documentation");
      hasTitle <- (regexpr("(@|[\\])title", bfr) != -1);
      if (!hasTitle)
        line <- paste(line, "\\title{", name, "}\n", sep="");
      objectName <<- value;
      usage <<- NULL;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocAbout <- function(bfr) {
      # An @RdocAbout does *not* have a value
      # bfr <- getTagValue(bfr);
      name <- ". About this package";
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "}\n", sep="");
      line <- paste(line, "\\alias{", name, "}\n", sep="");
      addKeyword("documentation");
      line <- paste(line, "\n", sep="");
      hasTitle <- (regexpr("(@|[\\])title", bfr) != -1);
      if (!hasTitle)
        line <- paste(line, "\\title{About this package}\n", sep="");
      objectName <<- name;
      usage <<- NULL;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagRdocPackage <- function(bfr) {
      # An @RdocPackage takes the package name as an argument
      bfr <- getTagValue(bfr);
      name <- attr(bfr, "value");
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "-package}\n", sep="");
      line <- paste(line, "\\alias{", name, "-package}\n", sep="");
      line <- paste(line, "\\alias{", name, "}\n", sep="");
      line <- paste(line, "\\docType{package}\n", sep="");
      addKeyword("package");
      line <- paste(line, "\n", sep="");
      hasTitle <- (regexpr("(@|[\\])title", bfr) != -1);
      if (!hasTitle)
        line <- paste(line, "\\title{Package ", name, "}\n", sep="");
      name <<- paste(name, "-package", sep="");
      objectName <<- name;
      usage <<- NULL;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagName <- function(bfr) {
      bfr <- getTagValue(bfr);
      if (!is.null(name)) {
        warning("Tag ignored: @RdocDefault is not needed if @RdocClass is specified.");
        return(bfr);
      }
      value <- attr(bfr, "value");
      name <- value;
      name <<- name <- escapeName(name);
      line <- paste("\\name{", name, "}\n", sep="");
      line <- paste(line, "\\alias{", name, "}", sep="");
      objectName <<- value;
      usage <<- Rdoc$getUsage(method=value);
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagAlias <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      alias <- value;
      alias <- escapeAlias(alias);
      line <- paste("\\alias{", alias, "}", sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagAliasMethod <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      method <- value;
      method <- escapeAlias(method);
      line <- paste("\\alias{", class, ".", method, "}\n", sep="");
      line <- paste(line, "\\alias{", method, ".", class, "}\n", sep="");
      line <- paste(line, "\\alias{", method, ",", class, "-method}\n", sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagAliasUndocumented <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      alias <- unlist(tools::undoc(package=value));
      alias <- sapply(alias, FUN=escapeAlias);
      lines <- paste("\\alias{", alias, "}", sep="");
      lines <- paste(lines, collapse="\n");
      rd <<- paste(rd, lines, sep="");
      bfr;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagUsage <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      if (nchar(value) == 0L) {
        # Usage is inferred from the @RdocNnn tag.
      } else if (nchar(value) > 0L) {
        # Usage is inferred from the @usage <spec> tag.
        parts <- strsplit(value, split=",", fixed=TRUE)[[1]];
        nparts <- length(parts);
        method <- parts[1L];
        if (nparts == 1L) {
          usage <- Rdoc$getUsage(method=method);
        } else if (nparts == 2L) {
          class <- parts[2L];
          usage <- Rdoc$getUsage(method=method, class=class);
        }
      }
      usage <- paste(usage, collapse="\n");
      line <- usage;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    tagSynopsis <- function(bfr) {
      usage <- c("", usage, "");
      usage <- paste(usage, collapse="\n");
      line <- paste("\\usage{", usage, "}", sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagKeyword <- function(bfr) {
      bfr <- getTagValue(bfr);
      keyword <- attr(bfr, "value");
      addKeyword(keyword);
      bfr;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagInclude <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      if (!file.exists(value)) {
        throw(RdocException("File to be included not found: ", value, source=sourcefile));
      } else {
        include <- readLines(value, warn=FALSE);
        include <- paste(include, collapse="\n");
        include <- gsub("\\%", "\\\\%", include);
        line <- paste(include, "\n", sep="");
        rd <<- paste(rd, line, sep="");
      }
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagExamples <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      if (!file.exists(value)) {
        throw(RdocException("File containing examples to be included not found: ", value, source=sourcefile));
      } else {
        include <- readLines(value, warn=FALSE);
        include <- paste(include, collapse="\n");
        include <- gsub("\\%", "\\\\%", include);
        line <- paste("\\examples{\n", include, "\n}", sep="");
      }
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagVisibility <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      visibility <<- value;
      if (!is.element(visibility, c("private", "protected", "public")))
        throw(RdocException("Unknown type of visibility: ", value, source=sourcefile));
      if (visibility == "private")
        addKeyword("internal");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagGet <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      if (exists(value)) {
        line <- as.character(get(value));
      } else if (!is.null(getOption(value))) {
        line <- as.character(getOption(value));
      } else {
        throw(RdocException("R variable does not exist: ", value, source=sourcefile));
      }
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagSet <- function(bfr, envir=parent.frame()) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      keyValue <- strsplit(value, "=")[[1]];
      assign(keyValue[1], value=keyValue[2], envir=envir);
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagAuthor <- function(bfr) {
      bfrT <- getTagValue(bfr);
      value <- attr(bfrT, "value");
      value <- as.character(value);
      hasValue <- (nchar(value) > 0L);
      hasValue <- hasValue && (regexpr("^[ \t]*[\n\r]", value) == -1L);

      # Does the @author tag has a value?
      if (hasValue) {
        # Non-empty @author tag with value, e.g. '@author "HB"'

        value <- gsub("^[ \t]*['\"]?", "", value);
        value <- gsub("['\"]?[ \t]*$", "", value);

        # (i) All authors?
        if (value == "*") {
          pkgAuthors <<- authors <- getPackageAuthors();
        } else {
          # (ii) All initials?  An initial = 2-5 upper case letters
          tmp <- unlist(strsplit(value, split=",", fixed=TRUE))
          tmp <- gsub("^[ \t]*", "", tmp);
          tmp <- gsub("[ \t]*$", "", tmp);
          tmpU <- toupper(tmp);
          pattern <- sprintf("^[%s]{2,5}$", paste(base::LETTERS, collapse=""));
          allInitials <- all( (tmpU == tmp) & (regexpr(pattern, tmp) != -1L) );
          if (allInitials) {
            initials <- tmp;

            # Create all initials of the 'authors'
            pkgAuthors <<- authors <- getPackageAuthors();
            fullnames <- format(authors, include=c("given", "family"));
            known <- abbreviate(fullnames, minlength=2L);
            known <- toupper(known);

            # Check if the initials match
            idxs <- match(initials, known);
            unknown <- initials[is.na(idxs)];
            if (length(unknown) > 0L) {
              known <- paste(sQuote(known), sQuote(fullnames), sep="=");
              throw(RdocException("Rdoc 'author' tag specifies initials (", paste(sQuote(unknown), collapse=", "), ") that are not part of the known ones (", paste(known, collapse=", "), ")", source=sourcefile));
            }
            authors <- authors[idxs];
          } else {
            authors <- as.person(value);
          }
        }
        bfr <- bfrT;
      } else {
        # Empty @author tag, i.e. '@author'

        pkgAuthors <<- authors <- getPackageAuthors();
        # If there are creators of the package (which there should be),
        # use those as the default for an empty '@author' tag.
        isCreator <- sapply(authors, FUN=function(a) is.element("cre", a$role));
        if (any(isCreator)) {
          authors <- authors[isCreator];
        }

        if (exists("author", mode="character", envir=globalenv())) {
          if (!authorWarn) {
            author <- get("author", mode="character", envir=globalenv());
            warning("Detected an 'author' character variable in the global environment. Note that, since R.oo 1.13.0, the author for an (empty) Rdoc @author tag is inferred from the 'Authors@R' or 'Author' DESCRIPTION field and no longer take from the global 'author' variable: ", sQuote(author));
            authorWarn <<- TRUE;
          }
        }
      }

      authorsTag <- format(authors, include=c("given", "family"));
      authorsTag <- paste(authorsTag, collapse=", ");

      line <- paste("\\author{", authorsTag, "}", sep="");
      rd <<- paste(rd, line, sep="");

      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagHowToCite <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      package <- Package(value);
      howToCite <- getHowToCite(package, newline=NULL);
      if (is.null(howToCite)) {
        line <- "\\emph{No citation information available.}\n";
      } else {
        line <- gsub("\n", " ", howToCite);
        line <- gsub("[ ]+", " ", line);
        line <- lapply(line, FUN=strwrap, width=85L);
        line <- lapply(line, FUN=paste, collapse="\n");
        line <- unlist(line, use.names=FALSE);
        line <- paste(line, collapse="\n\n");
        line <- paste("\\preformatted{\n", line, "\n}\n", sep="");
        # Add the following line to fix a "R CMD check-bug" in LaTeX.
        # /HB 2004-03-10
        line <- paste(line, "\\emph{}\n", sep="");
      }
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagSeeclass <- function(bfr) {
      line <- paste("For more information see \\code{\\link{", class, "}}.", sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagSeemethod <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      line <- paste("\\code{\\link[", package, ":", Rdoc$createName(class, value), "]{*", value, "}()}", sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagSeeOtherMethods <- function(bfr) {
      clazz <<- Class$forName(class);
      method <- gsub(paste("[.]", class, "$", sep=""), "", objectName);

      classes <- list(getSuperclasses(clazz), getKnownSubclasses(clazz));
      for (k in 1:2) {
        methods <- paste(method, classes[[k]], sep=".");
        exists <- sapply(methods, FUN=exists, mode="function");
        classes[[k]] <- classes[[k]][exists];
        methods <- methods[exists];
        for (l in seq_along(methods)) {
          fcn <- get(methods[l], mode="function");
          modifiers <- attr(fcn, "modifiers");
          isPrivate <- is.element("private", modifiers);
          isDeprecated <- is.element("deprecated", modifiers);
          if (isPrivate || isDeprecated)
            methods[l] <- NA;
        }
        classes[[k]] <- classes[[k]][!is.na(methods)];
      }

      for (k in 1:2) {
        classes[[k]] <- paste("\\code{\\link[", package, ":", Rdoc$createName(classes[[k]], method), "]{", classes[[k]], "}}", sep="");
      }

      line <- "\\cr\\bold{Implementations of this method in}";
      line <- paste(line, " i) \\bold{superclasses:} ");
      line <- paste(line, paste(classes[[1]], collapse=", "), sep="");
      line <- paste(line, ", and ii) \\bold{subclasses:} ");
      line <- paste(line, paste(classes[[2]], collapse=", "), ".", sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagSee <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      pkgObject <- strsplit(value, "::", value)[[1]];
      fcn <- "";
      if (length(pkgObject) == 1) {
        pkg <- NULL;
        obj <- pkgObject[1];
        if (exists(obj, mode="function")) {
          expr <- substitute(inherits(fcn, "Class"), list(fcn=as.name(obj)));
          if (!eval(expr))
            fcn <- "()";
        }
      } else {
        pkg <- pkgObject[1];
        obj <- pkgObject[2];
        .require <- require;  # To please R CMD check
        if (.require(package=pkg, character.only=TRUE)) {
          pos <- which(paste("package:", "base", sep="") == search());
          if (exists(obj, where=pos, mode="function", inherits=FALSE))
            fcn <- "()";
        }
        pkg <- paste("[", pkg, "]", sep="");
      }
      line <- paste("\\code{\\link", pkg, "{", obj, "}}", fcn, sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagClasshierarchy <- function(bfr) {
      txt <- getRdHierarchy(clazz);
      subclasses <- getKnownSubclasses(clazz);

      # If possible, create links to Rd docs for each of the subclasses.
      links <- c();
      for (name in subclasses) {
        link <- name;
        if (exists(name, mode="function")) {
          cls <- get(name, mode="function");
          if (inherits(cls, "Class")) {
            pkg <- getPackage(cls);
            if (is.null(pkg))
              link <- paste("\\link{", link ,"}", sep="")
            else
              link <- paste("\\link[", pkg, "]{", link ,"}", sep="");
            if (isAbstract(cls))
              link <- paste("\\emph{", link, "}", sep="");
          }
        }

#        link <- paste("\\code{", link ,"}", sep="");
        links <- c(links, link);
      } # for (name in ...)
      subclasses <- paste(links, collapse=", ");

      txt <- paste(txt, "\\bold{Directly known subclasses:}\\cr\n", sep="");
      txt <- paste(txt, subclasses, sep="");
      txt <- paste(txt, "\\cr\n\n", sep="");
      txt <- paste(txt, getRdDeclaration(clazz), "\n", sep="");
      line <- txt;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagMethods <- function(bfr) {
      methods <- getRdMethods(clazz, visibility=visibility);
      line <- methods;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagAllMethods <- function(bfr) {
      bfr <- getTagValue(bfr);
      visibilities <- attr(bfr, "value");
      # Ad hoc patch for parser /060530.
      addEnd <- (identical(visibilities, "}"))
      visibilities <- gsub(" ", "", visibilities);
      visibilities <- unlist(strsplit(visibilities, split="|", fixed=TRUE));
      methods <- getRdMethods(clazz, visibilities=visibilities);
      line <- paste(methods, "\n\n", sep="");
      methods <- Rdoc$methodsInheritedFrom(clazz, visibility, showDeprecated=showDeprecated, sort=sort);
      line <- paste(line, methods, sep="");
      if (addEnd)
        line <- paste(line, "}", sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagTitle <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      title <<- as.character(value);

      # Check the title for correctness according to http://developer.r-project.org/Rds.html
      firstLetter <- substring(title, 1,1);
      if (firstLetter != toupper(firstLetter))
        throw(RdocException("Titles shoule be capitalized: ", title, source=sourcefile));
      if (regexpr("[.]$", title) != -1)
        throw(RdocException("Titles should not end with a period: ", title, source=sourcefile));
      if (regexpr("[^\\][\\][:letter:]", title) != -1)
        throw(RdocException("Titles should not contain markup: ", title, source=sourcefile));

      if (isDeprecated)
        title <<- paste("Deprecated: ", title, sep="");
      line <- paste("\\title{", title, "}", sep="");
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagEval <- function(bfr) {
      bfr <- getTagValue(bfr);
      value <- attr(bfr, "value");
      expr <- as.character(value);
      result <- eval(parse(text=expr));
      line <- result;
      rd <<- paste(rd, line, sep="");
      bfr;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tagInheritedmethods <- function(bfr) {
      methods <- Rdoc$methodsInheritedFrom(clazz, visibility, sort=sort);
      line <- methods;
      rd <<- paste(rd, line, sep="");
      bfr;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tags <- list(
      "RdocClass"         = tagRdocClass,
      "RdocGeneric"       = tagRdocGeneric,
      "RdocFunction"      = tagRdocFunction,
      "RdocDefault"       = tagRdocDefault,
      "RdocMethod"        = tagRdocMethod,
      "RdocObject"        = tagRdocObject,
      "RdocData"          = tagRdocData,
      "RdocDocumentation" = tagRdocDocumentation,
      "RdocAbout"         = tagRdocAbout,
      "RdocPackage"       = tagRdocPackage,
  #
      "classhierarchy"   = tagClasshierarchy,  # must be *before* "class".
      "synopsis"         = tagSynopsis,
      "usage"            = tagUsage,
      "keyword"          = tagKeyword,
  #
      "Class"            = tagRdocClass,
      "methodname"       = tagRdocMethod,
      "name"             = tagName,
  #
      "aliasmethod"      = tagAliasMethod,     # must be *before* "alias".
  #   "aliasundocumented" = tagAliasUndocumented, # not useful.
      "alias"            = tagAlias,
      "title"            = tagTitle,
  #
      "class"            = tagClass,
      "include"          = tagInclude,
      "inheritedmethods" = "tagInheritedmethods",
      "examples"         = tagExamples,
      "visibility"       = tagVisibility,
      "get"              = tagGet,
      "set"              = tagSet,
      "author"           = tagAuthor,
      "howtocite"        = tagHowToCite,
      "seeclass"         = tagSeeclass,
      "seemethod"        = tagSeemethod,
      "seeothermethods"  = tagSeeOtherMethods,
      "see"              = tagSee,
      "methods"          = tagMethods,
      "allmethods"       = tagAllMethods,
      "eval"             = tagEval
    );
    names <- names(tags);
    attr(tags, "beginsWith") <- paste("^@", names, sep="");


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # SETUP
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Set the locale
    Sys.setlocale(locale=locale);


    # Make a best guess what the package is that is created by querying
    # the DESCRIPTION file, which should be in the parent directory of
    # the getManPath() directory.
    pkgPath <- dirname(getManPath(this));
    pathname <- file.path(pkgPath, "DESCRIPTION");
    if (!file.exists(pathname)) {
      stop("Cannot infer package name. File not found: ", pathname);
    }
    pi <- read.dcf(file=pathname);
    package <- pi[,"Package", drop=TRUE];
    if (length(package) == 0L) {
      throw("Failed to infer package name.  No 'Package' was listed in ", sQuote(pathname), ".");
    }
    if (length(package) > 1L) {
      throw("Failed to infer package name.  More than one 'Package' were listed in ", sQuote(pathname), ": ", paste(sQuote(package), collapse=", "));
    }

    Rdoc$package <- package;

    class <- NULL;
    clazz <- NULL;
    Rdoc$source <- sourcefile <<- NULL;

    rds <- list();
    for (rdoc in rdocs) {
      # Remember the name of the source file in case of an error...
      Rdoc$source <- sourcefile <<- attr(rdoc, "sourcefile");

      title <- NULL;
      objectName <- NULL;
      isDeprecated <- FALSE;
      method <- NULL;
      name <- NULL;
      usage <- NULL;
      visibility <- NULL;

      # ==============================================================
      # 1. Replace all shorttags
      #
      #    This can be vectorized since all shorttags are read only
      #    and does not rely on any state of Rdoc etc.
      #    This should speed up the process.
      # ==============================================================
      for (kk in seq_along(shorttags)) {
        replace <- attr(shorttags, "replace")[kk];

        # (a) Replace all occurances at the beginning of the lines.
        pattern <- attr(shorttags, "beginsWith")[kk];
        rdoc <- gsub(pattern, replace, rdoc);

        # (b) Replace all other occurances.
        ready <- FALSE;
        while (!ready) {
          pattern <- attr(shorttags, "contains")[kk];
          pos <- regexpr(pattern, rdoc);
          idx <- (pos != -1L);
          if (any(idx)) {
            len <- attr(pos, "match.length")[idx];
            pos <- pos[idx];
            prefix <- substring(rdoc[idx], first=1L, last=pos);
            suffix <- substring(rdoc[idx], first=pos+len);
            rdoc[idx] <- paste(prefix, replace, suffix, sep="");
          } else {
            ready <- TRUE;
          }
        }
      } # for (kk ...)


      # ==============================================================
      # 2. Compile the remaining lines
      #
      #    This can *not* be vectorized since some tags may change the
      #    state of the Rdoc compiler or set a local variable, which
      #    will affect following tags.
      # ==============================================================

      # (a) Make on big string of the whole Rdoc comment.
      #     This will simplify handling of line breaks within a tag
      #     argument, e.g. when @title "Bla bla\nbla" exists.
      rdoc <- paste(rdoc, collapse="\n");

      # (b) Get the first occurance of a tag.
      #     A tag begins with a '@', but can not have a letter, a
      #     number or a '.' in front since then it could be a code
      #     fragment refering to a S4 slot.
      patternL <- "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
      pattern <- sprintf("[^%s0-9.]@", patternL);
      pattern2 <- sprintf("[^%s0-9]@", patternL);

      rd <- "";
      while (TRUE) {
        # (i) It can be that the tag start at position 1...
        pos <- regexpr("^@", rdoc);

        # (ii) ...otherwise...
        if (pos == -1L) {
          pos <- regexpr(pattern, rdoc);
          if (pos == -1L)
            break;  # No more tags. Done!
          # Save everything before the tag...
          rd <- paste(rd, substring(rdoc, first=1L, last=pos), sep="");
          # ...and make the tag be at the first position.
          rdoc <- substring(rdoc, first=pos+1L);
        }

        if (debug)
          str(rdoc);

        # (iii) Identify the tag
        tagName <- NULL;
        tagFunction <- NULL;
        for (kk in seq_along(tags)) {
          patternKK <- attr(tags, "beginsWith")[kk];
          pos <- regexpr(patternKK, rdoc);
          if (pos != -1L) {
            len <- attr(pos, "match.length");
            tagIdx <- kk;
            tagName <- names(tags)[kk];
            if (debug)
              cat(paste("Found tag: ", tagName, "\n", sep=""));
            tagFunction <- tags[[kk]];
            break;
          }
        } # for (kk ...)

        if (!is.null(tagFunction)) {
          # Shift the Rdoc buffer
          rdoc <- substring(rdoc, first=len+1L);
          if (is.function(tagFunction)) {
            # Evaluate the tag function in the current environment
            # so all variables can be shared between tags.
            # All tag functions must return the resulting buffer!
            expr <- substitute(tagFunction(rdoc),
                      list(tagFunction=tagFunction, rdoc=rdoc));
            rdoc <- eval(expr);
          }
        } else {
          pos <- regexpr(pattern2, substring(rdoc, first=2L));
          tagName <- substring(rdoc, first=1L, last=pos);
          msg <- sprintf("Unknown tag not processed in '%s': '%s'", sourcefile, tagName);
          warning(msg);
          rd <- paste(rd, substring(rdoc, first=1L, last=1L), sep="");
          rdoc <- substring(rdoc, first=2L);
        }

        if (isDeprecated && !showDeprecated)
          break;
      } # while(TRUE), i.e. get first tag...

      if (showDeprecated || !isDeprecated) {
        # Do not forget to add the rest!
        rd <- paste(rd, rdoc, sep="");
        rdoc <- NULL;


        # Append all keywords at the end
        rd <- paste(rd, getRdKeywords(), sep="\n");

        # Remove all empty lines
        rd <- gsub("[ \t]\n", "\n", rd);
        rd <- gsub("[ \t]\r", "\r", rd);

        if (is.null(name)) {
          # @RdocClass, @RdocDefault and/or @RdocMethod was not given. Search for classical \name{...}
          search <- regexpr("\\name\\{[^\\}]*\\}", rd);
          if (search == -1L) {
            throw(RdocException("The resulting Rd text does not have a \\name{} tag: ", substring(rd, first=1L, last=40L), source=sourcefile));
          }
          name <- substring(rd, first=search+5L, last=search+attr(search, "match.length")-2);
          search <- regexpr("\\name\\{[^\\}]*\\}", substring(rd, first=search+1L));
          if (search != -1L)
            throw(RdocException("The resulting Rd text has more than one \\name{} tag.", source=sourcefile));
        }

        visibility <- "public";
        if (is.null(visibility)) {
          if (is.null(objectName)) {
          } else if (!exists(objectName)) {
            # If no object was found, assume that it is a Rdoc comment for
            # a non-object, i.e. a concept or similar.
          } else {
            object <- get(objectName);
            modifiers <- attr(object, "modifiers");
            if (is.element("private", modifiers)) {
              visibility <- "private";
            } else if (is.element("protected", modifiers)) {
              visibility <- "protected";
            }
          }
        }

        attr(rd, "visibility") <- as.character(visibility);
        attr(rd, "isDeprecated") <- isDeprecated;
        attr(rd, "name") <- as.character(name);
        attr(rd, "sourcefile") <- sourcefile;

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Check Rd code?
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if (check) {
          if (compareVersion(as.character(getRversion()), "2.10.0") >= 0) {
            # R v2.10.0 and newer
            tryCatch({
              con <- textConnection(rd);
              rdParse <- tools::parse_Rd(file=con);
            }, warning = function(w) {
              filename <- sprintf("%s.Rd.ERROR", attr(rd, "name"));
              cat(rd, sep="\n", file=filename);
              throw(RdocException(sprintf("Syntax error in generated Rd code (see '%s') for Rdoc comment '%s' (in '%s') was detected by tools:parse_Rd(): %s", filename, attr(rd, "name"), attr(rd, "sourcefile"), as.character(w))));
            }, finally = {
              close(con);
              con <- NULL;
            })
          } else {
            # R v2.9.2 and before
            tryCatch({
              ns <- getNamespace("tools");
              tools_Rd_parse <- get("Rd_parse", mode="function", envir=ns);
              rdParse <- tools_Rd_parse(text=rd);
              if (length(rdParse$rest) > 0) {
                throw(RdocException("Unknown top-level text in generated Rd code for Rdoc comment '", attr(rd, "name"), "' (in '", attr(rd, "sourcefile"), "') (typically due to too many or a missing bracket): ", paste(rdParse$rest, collapse=", ", sep="")));
              }
            }, error = function(e) {
              throw(RdocException("Syntax error in generated Rd code for Rdoc comment '", attr(rd, "name"), "' (in '", attr(rd, "sourcefile"), "') was detected by tools:Rd_parse(): ", as.character(e)));
            })
          }
        } # if (check)

        rds <- c(rds, list(rd));
      } else {
        warning(paste("No Rd file for '", objectName, "' was generated since it was declared deprecated.", sep=""));
      } # if (!isDeprecated)
    } # for (rdoc in rdocs)

    rds;
  } # compileRdoc()


  # A local version of the sourceTo() in R.io.
  sourceTo <- function(..., local=TRUE, envir=parent.frame()) {
    # Wrap up the arguments to source
    args <- list(...);

    if (!is.element("file", names(args)))
      names(args)[1] <- "file";

    # Override any 'local' argument
    args$local <- local;

    # Create a call expression to source(..., local=local)
    expr <- substitute({
      do.call(source, args)
    }, list(args=args));

    # Call source()
    eval(expr, envir=envir);
  } # sourceTo()
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Assert that the manPath exists.
  createManPath(this);

  filename <- as.character(filename);
  if (length(filename) == 1L && !file.exists(filename)) {
    if (regexpr("\\*", filename) != -1L || regexpr("\\?", filename) != -1L) {
      # Filename pattern
      pattern <- filename;
      # List all files
      filename <- list.files();
      # Match to filename pattern
      filename <- grep(paste(pattern, "$", sep=""), filename, value=TRUE);
      # Keep only files
      filename <- filename[file_test("-f", filename)];
      # Keep only newer files?
    }
  }

  # Load the source code in case it contains new stuff.
  if (source) {
    lapply(filename, FUN=source);
  }

  for (file in filename) {
    if (verbose)
      cat("Generating Rd files from source file ", file, "...", sep="");

    rdocs <- extractRdocs(file, verbose=verbose, debug=debug);
    rd <- compileRdoc(rdocs, showDeprecated=showDeprecated, verbose=verbose, debug=debug);
    writeRd(rd, path=destPath, addTimestamp=addTimestamp, verbose=verbose, debug=debug);
    if (verbose)
      cat("\n");
  }

  if (check) {
    check(this, manPath=destPath, verbose=verbose);
  }
}, static=TRUE) # compile()







###########################################################################/**
# @RdocMethod hierarchy
#
# @title "Gets the class hierarchy"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("hierarchy", "Rdoc", function(this, class, ...) {
  package <- getPackage(class);
  s <- paste("Package: ", package, "\\cr\n");
  what <- if (inherits(class, "Class")) "Class" else "Interface";
  s <- paste(s, "\\bold{", what, " ", getName(class), "}\\cr\n\n", sep="");
  indent <- "";
  for (extend in rev(getSuperclasses(class))) {
    link <- sapply(extend, FUN=function(name) {
#      isAbstract <- FALSE;
      link <- name;
      if (exists(name, mode="function")) {
        cls <- get(name, mode="function");
        if (inherits(cls, "Class")) {
          pkg <- getPackage(cls);
          if (is.null(pkg))
            link <- paste("\\link{", link ,"}", sep="")
          else
            link <- paste("\\link[", pkg, "]{", link ,"}", sep="");
#          if (isAbstract(cls)) {
#            link <- paste("\\emph{", link, "}", sep="");
#            isAbstract <- TRUE;
#          }
        }
      }
      paste("\\code{", link ,"}", sep="");
    });
    if (indent == "") {
      s <- paste(s, link, "\\cr\n", sep="");
      indent <- "~~";
    } else {
      s <- paste(s, "\\code{", indent, "+--}", link, "\\cr\n", sep="");
      indent <- paste(indent, "~~~~~", sep="");
    }
    s <- paste(s, "\\code{", indent, "|}\\cr\n", sep="");
  }
  link <- paste("\\code{", getName(class), "}", sep="");
  if (isAbstract(class))
    link <- paste("\\emph{", link, "}", sep="");
  s <- paste(s, "\\code{", indent, "+--}", link, "\\cr\n\n", sep="");
  s;
}, private=TRUE, static=TRUE);






###########################################################################/**
# @RdocMethod declaration
#
# @title "Gets the class declaration"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("declaration", "Rdoc", function(this, class, ...) {
  s <- "public"; # visibility(class);
  if (isAbstract(class))
    s <- paste(s, "abstract");
  if (isStatic(class))
    s <- paste(s, "static");
  if (inherits(class, "Class"))
    s <- paste(s, "class")
  else
    throw(getName(class), " is neither a class nor an interface.");

  s <- paste(s, " \\bold{", getName(class), "}\\cr\n", sep="");
  links <- getSuperclasses(class);

  if (length(links) > 0) {
    name <- links[1];
    link <- name;
    if (exists(name, mode="function")) {
      cls <- get(name, mode="function");
      if (inherits(cls, "Class")) {
        pkg <- getPackage(cls);
        if (is.null(pkg))
          link <- paste("\\link{", link ,"}", sep="")
        else
          link <- paste("\\link[", pkg, "]{", link ,"}", sep="");
        if (isAbstract(cls))
          link <- paste("\\emph{", link, "}", sep="");
      }
    }
    paste("\\code{", link ,"}", sep="");
    s <- paste(s, "extends ", link, "\\cr\n", sep="");
  }
  s;
}, private=TRUE, static=TRUE);





###########################################################################/**
# @RdocMethod methodsInheritedFrom
#
# @title "Gets all methods inherited from a class in Rd format"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("methodsInheritedFrom", "Rdoc", function(this, class, visibility=c("public", "protected", "private"), showDeprecated=FALSE, inheritedFrom=NULL, sort=TRUE, trial=FALSE, ...) {
  s <- "";

  private <- is.element("private", visibility);

  # Classes
  for (extend in getSuperclasses(class)) {
    # Try to find a Class object with this name.
    clazz <- NULL;
    tryCatch({
      clazz <- Class$forName(extend);
    }, error = function(ex) {})

    if (is.null(clazz)) {
      # Use methods() to find methods
      methods <- methods(class=extend);
      pattern <- paste("[.]", extend, "$", sep="");
      methods <- gsub(pattern, "", methods);
    } else {
      # Get all methods of this Class
      methods <- getMethods(clazz, private=private, deprecated=showDeprecated);
      methods <- methods[[extend]];
      methods <- names(methods);
    }
    if (length(methods) > 0) {
      methods <- paste(methods, collapse=", ");
      s <- paste(s, sprintf("\\bold{Methods inherited from %s}:\\cr\n", extend));
      s <- paste(s, methods, "\n\n", sep="");
    }
  }

  s;
}, private=TRUE, static=TRUE);



setMethodS3("getObject", "Rdoc", function(static, name, mode="any", package=static$package, ...) {
  # Search for object in:
  #  (i) the package names iff found, then
  # (ii) then the search path.

  # Try to infer the package's namespace.
  ns <- NULL;
  if (is.character(package)) {
    tryCatch({
      ns <- getNamespace(package);
    }, error = function(ex) {});
  }

  if (is.environment(ns) && exists(name, mode=mode, envir=ns)) {
    obj <- get(name, mode=mode, envir=ns);
  } else if (exists(name, mode=mode)) {
    obj <- get(name, mode=mode);
  } else {
    throw("Failed to locate object of mode \"", mode, "\": ", name);
  }

  obj;
}, private=TRUE, static=TRUE);



###########################################################################/**
# @RdocMethod getUsage
#
# @title "Gets the usage of a method"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{method}{A method name (@character string).}
#   \item{class}{An optional class name (@character string).}
#   \item{wrap}{An @integer specifying the maximum number of characters per line.  Longer lines will be wrapped with newlines.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getUsage", "Rdoc", function(static, method, class=NULL, wrap=90L, indent=2L, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  buildUsage <- function(method, class=NULL, args, valueArg=NULL, wrap=90L, head=NULL, ...) {
    # Argument 'args':
    stopifnot(is.list(args));

    indentStr <- paste(rep(" ", times=indent), collapse="");

    correction <- 0L;
    if (length(head) == 0L) {
      if (is.null(class)) {
        head <- method;
        # Escape '%*%' to '\%*\%'
        head <- gsub("%", "\\%", head, fixed=TRUE);
        # Quote any method name containing '%'
        if (regexpr("%", head, fixed=TRUE) != -1L)
          head <- sprintf("`%s`", head);
      } else {
        # The effective length of this in the help manual is nchar(method).
        head <- sprintf("\\method{%s}{%s}", method, class);
        correction <- nchar(head) - nchar(method);
      }
    }
    correction0 <- correction;

    lines <- NULL;
    line <- paste(head, "(", sep="");
    if (length(args) == 0L) {
      line <- paste(line, ") ", sep="");
    }

    while (length(args) > 0L) {
      subargs <- args[[1L]];
      nsubargs <-length(subargs);

      # Try to keep <key>=<value> together
      if (nsubargs >= 3L) {
        # If <key>=<value> fit on a line, then keep the together...
        if (sum(nchar(subargs[1:3L])) <= wrap) {
          subargs[3L] <- paste(subargs[1:3], collapse="");
          subargs <- subargs[-(1:2)];
        } else if (sum(nchar(subargs[1:2L])) <= wrap) {
        # ...otherwise, at least keep <key>= together, iff possible.
          subargs[2L] <- paste(subargs[1:2], collapse="");
          subargs <- subargs[-1L];
        }
        nsubargs <-length(subargs);
      }

      # Remaining arguments
      args <- args[-1L];
      nargs <- length(args);
      suffix <- if (nargs > 0L) ", " else ") ";

      # For each subargument
      for (kk in seq_len(nsubargs)) {
        subarg <- subargs[kk];
##        str(list(kk=kk, subarg=subarg));

        if (kk == nsubargs) {
          subarg <- paste(subarg, suffix, sep="");
        } else {
          subarg <- paste(subarg, " ", sep="");
        }
        len <- nchar(subarg);

        # Does argument fit on the same line?
        if (nchar(line) - correction + len <= wrap) {
          line <- paste(line, subarg, sep="");
        } else {
          lines <- c(lines, line);
          line <- paste(indentStr, subarg, sep="");
          correction <- 0L;
        }
      } # for (kk ...)
    } # while (length(args) > 0L)

    # Append a value assignment, i.e. "... <- value"?
    if (!is.null(valueArg)) {
      arg <- paste("<- ", valueArg, sep="");
      # Does it fit on the same line?
      if (nchar(line) - correction + nchar(arg) <= wrap) {
        line <- paste(line, arg, sep="");
      } else {
        lines <- c(lines, line);
        line <- paste(indentStr, arg, sep="");
      }
    }
    lines <- c(lines, line);
    lines <- gsub("[ ]$", "", lines); # Trim trailing space

    # Sanity check
    lens <- nchar(lines);
    lens[1L] <- lens[1L] - correction0;
    stopifnot(all(lens <= wrap));

##    print(lines);

    lines;
  } # buildUsage()


  if (!is.null(class)) {
    fcnName <- paste(method, class, sep=".");
  } else {
    fcnName <- method;
  }

  fcn <- NULL;
  tryCatch({
    fcn <- getObject(static, name=fcnName, mode="function");
  }, error = function(ex) {
    cat("Failed...\n");
    print(ex);
    cat("Failed...done\n");
  });
  if (!is.function(fcn)) {
    throw(RdocException("Could not get usage. Function was not found: ", fcnName, "()", source=Rdoc$source));
  }

  isStatic <- is.element("static", attr(fcn, "modifiers"));
  isConstructor <- inherits(fcn, "Class");
  args <- Rdoc$argsToString(fcn, escapeRd=TRUE, collapse=FALSE);

  # Replacement methods are special
  isReplacement <- (regexpr("<-$", method) != -1L);
  if (isReplacement) {
    method <- gsub("<-$", "", method);
    nargs <- length(args);
    valueArg <- args[nargs];
    args <- args[-nargs];
  } else {
    valueArg <- NULL;
  }

  if (isConstructor) {
    usage <- buildUsage(method, args=args, valueArg=valueArg, wrap=wrap);
  } else if (isStatic) {
    # (a) The S3 method call
    lines <- buildUsage(method, class=class, args=args, valueArg=valueArg, wrap=wrap);
    usageM <- paste(lines, collapse="\n");

    # (b) The "static" method call, e.g. Class$forName(...)
    # Adjust line width ('wrap') to fit prefix '## ' as well.
    lines <- buildUsage(method, class=class, args=args[-1L], valueArg=valueArg, head=paste(class, method, sep="$"), wrap=wrap - 3L);
    lines <- paste("## ", lines, sep="");
    usageS <- paste(lines, collapse="\n");

    # (c) Combine
    usage <- c("## Static method (use this):",
               usageS,
               "",
               "## Don't use the below:",
               usageM);
  } else if (!is.null(class)) {
    usage <- buildUsage(method, class=class, args=args, valueArg=valueArg, wrap=wrap);
  } else {
    usage <- buildUsage(method, args=args, valueArg=valueArg, wrap=wrap);
  }

  usage;
}, private=TRUE, static=TRUE) # getUsage()




###########################################################################/**
# @RdocMethod getClassS4Usage
#
# @title "Gets the usage of a S4 class"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{class}{A class name (@character string).}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getClassS4Usage", "Rdoc", function(static, class, ...) {
  if (!inherits(class, "classRepresentation"))
    throw(InternalErrorException("Wrong usage."));

  name <- class@className;

  usage <- paste("new(", name, ")", sep="");

  hasConstructor <- exists(name, mode="function");
  if (hasConstructor) {
    constructor <- get(name, mode="function");
    args <- Rdoc$argsToString(constructor, collapse=TRUE);
    args <- paste(args, collapse=", ");
    constructorUsage <- paste(name, "(", args, ")", sep="");
    usage <- paste(usage, "\n", constructorUsage, sep="");
  }

  usage;
}, private=TRUE, static=TRUE);




###########################################################################/**
# @RdocMethod argsToString
#
# @title "Gets the arguments signature of a function"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{fcn}{A @function.}
#  \item{escapeRd}{If @TRUE, certain Rd markup symbols are escaped.}
#  \item{collapse}{If @TRUE, each argument is returned as a single string,
#   otherwise split up into a vector of strings as far as possible.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @list of @character strings.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("argsToString", "Rdoc", function(static, fcn, escapeRd=FALSE, collapse=TRUE, ...) {
  a <- args(fcn);

  # Nothing to do?
  if (is.null(a)) {
    return("[primitive function]");
  }

  # Sanity check
  if (typeof(a) != "closure") {
    throw("Expected closure but found something else: ", typeof(a));
  }

  args <- formals(a);
  argsNames <- names(args);

  res <- list();
  for (kk in seq_along(args)) {
    arg     <- args[kk];
    argName <- argsNames[kk];

    s <- argName;

    argDefault <- arg[[1L]];
    if (!missing(argDefault)) {
      argDefault <- deparse(argDefault, width.cutoff=20L);
      argDefault <- trim(argDefault);

      # Escape '%' (which is a comment in Rd format)?
      if (escapeRd) {
        argDefault <- gsub("\\%", "\\\\%", argDefault);
      }

      if (collapse) {
        argDefault <- paste(argDefault, collapse=" ");
        s <- paste(s, "=", argDefault, sep="", collapse="");
      } else {
        s <- c(s, "=", argDefault);
      }
    }

    res <- c(res, list(s));
  }

  res;
}, private=TRUE, static=TRUE) # argsToString()




###########################################################################/**
# @RdocMethod getRdTitle
#
# @title "Extracts the title string of a Rd file"
#
# \description{
#   @get "title" corresponding the the specified method of the specified class.
# }
#
# @synopsis
#
# \arguments{
#   \item{method}{The method to be search for.}
#   \item{class}{The class the method belongs to.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getRdTitle", "Rdoc", function(this, class, method, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tools_fetchRdDB <- get("fetchRdDB", mode="function", envir=getNamespace("tools"), inherits=FALSE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # e s c a p e N a m e ( )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  escapeName <- function(name) {
    name <- gsub("\\#", "POUND", name);
    name <- gsub("\\$", "DOLLAR", name);
    name <- gsub("\\%", "PERCENT", name);
    name <- gsub("\\&", "AND", name);
    name <- gsub("\\~", "TILDE", name);
    name <- gsub("\\_", "UNDERSCORE", name);
    name <- gsub("\\^", "POWER", name);
    name <- gsub("\\\\", "BACKSLASH", name);
    name <- gsub("\\{", "LCURLY", name);
    name <- gsub("\\}", "RCURLY", name);
    name <- gsub("<-", "< -", name);
    name;
  } # escapeName()

  title <- NULL;

  # Search for the file <class>.<method>.Rd in the man/ directory
  name <- createName.Rdoc(NULL, getName(class), method, escape=FALSE);
  name <- escapeName(name);

  rdName <- Rdoc$escapeRdFilename(name);
  rdFile <- paste(rdName, "Rd", sep=".");
  url <- file.path(getManPath(this), rdFile);
  if (file.exists(url)) {
    src <- paste(readLines(url, warn=FALSE), collapse="\n");

    # Search for \title{...} in the Rd source
    titlePos <- regexpr("\\title\\{[^\\}]*}", src);
    if (titlePos == -1) {
      warning(paste("Could not find a \\title{} definition in the Rd file for ", method, " in ", getName(class), ". Will search in loaded packages.", sep=""));
      "";
    } else {
      title <- trim(substring(src, first=titlePos+6, last=titlePos+attr(titlePos, "match.length")-2));
    }
  } else {
    warning(paste("The Rd file for ", method, " in ", getName(class), " could not be found. Will search in loaded packages.", sep=""));
  }

  if (is.null(title)) {
    methodName <- paste(method, ".", getName(class), sep="");
    packageName <- Rdoc$getPackageNameOf(methodName, mode="function");
    if (length(packageName) == 1) {
      if (compareVersion(as.character(getRversion()), "2.10.0") >= 0) {
        # R v2.10.0 and newer
        path <- system.file("help", package=packageName);
        filebase <- file.path(path, packageName);
        tryCatch({
          entry <- tools_fetchRdDB(filebase, key=methodName);
          tags <- lapply(entry, FUN=attr, "Rd_tag");
          idx <- which(tags == "\\title");
          if (length(idx) > 1) {
            idx <- idx[1];
          }
          if (length(idx) == 1) {
            entry <- entry[[idx]];
            entry <- entry[[1]];
            title <- entry[1];
          }
        }, error = function(ex) {
          warning(as.character(ex));
        });
      } else {
        package <- Package(packageName);
        tryCatch({
          contents <- getContents(package);
          pos <- which(contents[,"Entry"] == name);
          if (length(pos) == 0) {
            warning(paste("Reverted to the CONTENTS file of package '", packageName, "', but found not matching entry: ", name, sep=""));
          } else if (length(pos) > 2) {
            warning(paste("Found more than one matching entry in the CONTENTS file of package '", packageName, "'. Using the first one only: ", name, sep=""));
            pos <- pos[1];
          }
          if (length(pos) != 0) {
            title <- as.character(contents[pos, "Description"]);
            attr(title, "package") <- packageName;
          }
        }, error=function(ex) {
          warning(as.character(ex));
        })
      }
    }
  }

  title;
}, private=TRUE, static=TRUE);




###########################################################################/**
# @RdocMethod getPackageNameOf
#
# @title "Gets the package of a method or an object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{objectName}{An method or object name (@character string).}
#  \item{mode}{Mode of object (@character string).}
#  \item{unique}{If @TRUE, only the first matching package is returned if
#    more than one is found.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getPackageNameOf", "Rdoc", function(static, objectName, mode="any", unique=TRUE, ...) {
  # Search all namespaces that are *attached*
  pkgs <- grep("^package:", search(), value=TRUE)
  pkgs <- gsub("^package:", "", pkgs)
  found <- sapply(pkgs, FUN=function(pkg) {
    exists(objectName, mode=mode, envir=asNamespace(pkg))
  })
  package <- names(found)[found]
  if (length(package) == 1L) return(package)
  if (length(package) > 1L && unique) {
    warning("Found more than one occurance of '", objectName, "' among the attached namespaces. Will only return the first one: ", paste(sQuote(package), collapse=", "))
    return(package[1L])
  }

  # If not found, then search any other namespace *loaded*
  pkgs <- setdiff(loadedNamespaces(), pkgs)
  found <- sapply(pkgs, FUN=function(pkg) {
    exists(objectName, mode=mode, envir=asNamespace(pkg))
  })
  package <- names(found)[found]
  if (length(package) == 1L) return(package)
  if (length(package) > 1L && unique) {
    warning("Found more than one occurance of '", objectName, "' among the loaded namespaces. Will only return the first one: ", paste(sQuote(package), collapse=", "))
    return(package[1L])
  }

  character(0L)
}, private=TRUE, static=TRUE)



###########################################################################/**
# @RdocMethod check
#
# @title "Checks the compiled Rd files"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{manPath}{The path to the Rd files (@character string).}
#  \item{verbose}{If @TRUE, extra information is outputted.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a printable object, which, if non-empty will show the errors.
# }
#
# \details{
#  Internally the \code{tools} package is used.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("check", "Rdoc", function(this, manPath=getManPath(this), verbose=FALSE, ...) {
  # file paths with trailing '/' are not recognized! /HB 2004-10-13
  manPath <- gsub("/$", "", manPath);
  if (verbose)
    cat("Checking Rd files in '", manPath, "'...\n", sep="");

  if (compareVersion(as.character(getRversion()), "2.10.0") >= 0) {
    # For R (>= 2.10.0)
    pathnames <- list.files(pattern="[.]Rd$", path=manPath, full.names=TRUE);
    res <- NULL;
    for (kk in seq_along(pathnames)) {
      pathname <- pathnames[kk];
      res <- tools::checkRd(pathname);
    }
  } else {
    # For R (< 2.10.0)
    tools_check_Rd_files_in_man_dir <- get("check_Rd_files_in_man_dir", mode="function", envir=getNamespace("tools"), inherits=FALSE);
    res <- tools_check_Rd_files_in_man_dir(manPath);
    if (length(res$files_with_surely_bad_Rd) > 0) {
      throw("Syntax error in Rd file(s): ",
                          paste(res$files_with_surely_bad_Rd, collapse=", "));
    }

    if (length(res$files_with_likely_bad_Rd) > 0) {
      print(res$files_with_likely_bad_Rd);
      throw("Syntax error in Rd file(s): ",
                          paste(res$files_with_surely_bad_Rd, collapse=", "));
    }
  }

  if (verbose)
    cat("Checking Rd files in '", manPath, "'...done\n", sep="");
  res;
})



###########################################################################/**
# @RdocMethod isVisible
#
# @title "Checks if a member is visible given its modifiers"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{modifiers}{A @character string of modifiers.}
#  \item{visibilities}{A @character string of visibility flags.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns @TRUE if the modifiers are equal or higher than the visibility
#  flags, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("isVisible", "Rdoc", function(static, modifiers, visibilities, ...) {
  if (is.element("deprecated", modifiers) && !is.element("deprecated", visibilities))
    return(FALSE);

  if (is.element("trial", modifiers) && !is.element("trial", visibilities))
    return(FALSE);

  levels <- c("private", "protected", "public");
  modifiers <- intersect(modifiers, levels);
  if (length(modifiers) == 0)
    return(TRUE);

  visibilities <- intersect(visibilities, levels);
  if (length(visibilities) == 0)
    return(TRUE);

  modifiers <- factor(modifiers, levels=levels);
  visibilities <- factor(visibilities, levels=levels);

  any(as.integer(visibilities) <= as.integer(modifiers));
}, static=TRUE, protected=TRUE) # isVisible()



#########################################################################
# HISTORY:
# 2014-10-18
# o Now Rdoc$getPackageNameOf() also finds non-exported objects, and
#   it first searches all attached and then all loaded namespaces.
# 2014-04-26
# o Now Rdoc$getRdUsage() escapes '%*%' to '\%*\%'.
# 2013-10-07
# o Now Rdoc tag @howToCite does a better job when there are
#   multiple citations in CITATION.
# 2013-06-27
# o Added trial version of Rdoc tag @usage.
# 2013-05-30
# o Now Rdoc$compile() infer the package name from the DESCRIPTION
#   file (instead from the package directory name).
# 2013-05-19
# o Now Rdoc$getUsage() inserts line breaks so that any usage line
#   is at most 90 characters long.
# 2013-04-08
# o Now the @RdocData Rdoc tag also adds an \docType{data} Rd tag.
# 2013-04-04
# o BUG FIX: In R.oo v1.13.1 a bug was introduced causing @allmethods
#   to ignore the optional following tag value.
# 2013-04-03
# o Now Rdoc$compile() always outputs Rd files with '\n' line breaks
#   regardless of system.
# o Now Rdoc$compile() uses locale "C" (default) to assert that
#   the same Rd files are generated regardless of system settings.
# o Now Rdoc$compile() no longer gives a warning if argument 'filename'
#   is a vector.  It should also handle multiple filenames and
#   source=TRUE.
# 2013-03-25
# o BUG FIX: getTagValue() could read/consume a following '}', which
#   should not be considered a tag value.
# 2013-03-08
# o Added support for @author "John Doe" as well as @author "JD" where
#   the initials are then inferred from the package's DESCRIPTION file.
# 2012-12-28
# o Replaced all data.class(obj) with class(obj)[1].
# 2012-06-11
# o BUG FIX/GENERALIZATION: Rdoc$getKeywords() now uses system
#   environment variable R_DOC_DIR for locating the internal
#   KEYWORDS.db.  Thanks to Charles Hogg at NIST for suggesting this.
# 2012-04-17
# o Now Rdoc$getUsage() searches also the package namespace for
#   the function definition.  This is done, before searching the
#   global environment.
# o Added private static getObject().
# 2011-07-27
# o CLEANUP: Replaced a (!is.null(foo) > 0) with (!is.null(foo)).
# 2010-09-22
# o Now Rdoc lines are allowed to start with double ('##') or tripple
#   ('###') comment characters in addition to single ('#') ones.
# 2010-06-04
# o Now argument 'addTimestamp' of Rdoc$compile() default to FALSE.
#   This way the generate Rd file will remain identical unless there
#   are real Rdoc/code changes.  Not adding timestamps is better when
#   working with a version control systems.
# 2010-06-01
# o BUG FIX: If there are no Rd files, then check() of Rdoc would
#   throw the error "object 'res' not found".
# 2009-10-26
# o BUG FIX: Rdoc$compile() did not work with R v2.10.0 and newer.
# 2008-08-11
# o Replace all 'a %in% b' with is.element(a,b) due to weird bug, cf.
#   my R-devel post 'Argument "nomatch" matched by multiple actual
#   arguments ... %in% -> match?!?' on March 6, 2008.
# 2008-07-30
# o Added shorttags: 'raw'.
# 2007-09-17
# o Added 'warn=FALSE' to all readLines().
# 2007-06-09
# o Removed (incorrect) argument name 'list' from all substitute() calls.
# 2007-01-06
# o Now getMethodsInheritedFrom() recognizes visiblity private.
# 2006-09-12
# o Prepared the Rdoc compiler to no longer generating the \synopsis{}
#   statement for the @synopsis tag, which was the case for static
#   methods.  I got an early note from Kurt Hornik saying all R base and
#   recommended package will have \synopsis{} removed in favor of the
#   \usage{} statement by the release of R v2.4.0.  The Rd tag will be
#   deprecated around R v3.0.0 or so.
# 2006-05-29
# o Added protected static method isVisible().
# o Added argument to @allmethods to specify visibility etc.
# 2006-04-10
# o Replace the generated Rd comment header to exclude my name and the
#   date of the R.oo package, e.g. "2001-2006".  This was done to help
#   any diff.
# o BUG FIX: Rdoc$compile() did no longer write the name of the source
#   file in the header.
# o BUG FIX: The code for formatting replacement methods contained a bug
#   that generated an error.
# o Replaced "\t" with "\\t" in warning "Could not find a \title{}".
# o Added argument 'addTimestamp=TRUE' to Rdoc$compile().  This makes it
#   possible to turn of the timestamp, because timestamps makes diff,
#   say the one in Subversion, think there is a real different.
# 2006-04-03
# o Now replacement methods are recognized and treated specially, e.g.
#   \method{[}{MyClass}(i, j) <- value.
# 2006-03-14
# o BUG FIX: Compiling Rdoc comments with invalid keyword tags would
#   generate an internal error.  Same for invalid visibility tags etc.
# 2006-02-18
# o Now createName() escapes name by default.
# 2006-02-03
# o Added new tag RdocPackage.  This is to replace the ".About ..."
#   page, especially since <pkg>-package now is at the very top of
#   the package's index page.
# o An Rd filename must begin with a letter or a digit.  In R v2.2.1
#   and before this has not been a problem, but R CMD check in R v2.3.0
#   complains (with another error actually).
# 2006-01-10
# o Updated the date string in the generated Rd headers.
# 2006-01-06
# o Added Rd links to classes listed under "Directly known subclasses:".
# 2005-06-17
# o BUG FIX: Used invalid regular expression '\\name{[^\\}]*}' instead
#   of '\\name\\{[^\\}]*\\}'. Why? The correct string is literally
#   '\name\{[^\}]\}' (no escape codes). The '\{' and '\}' (no escape)
#   are to tell regexpr that it should match '{' and '}'; '{' and '}'
#   are used for different purposes. Simple, ehe. Thanks Lorenz Wernisch,
#   School of Crystallography, University of London of reporting this.
# 2005-06-08
# o BUG FIX: getRdKeywords() gave an error if no keywords are available.
# 2005-06-03
# o Added internal addKeyword() to Rdoc$compile().
# o Now an RdocMethod tag will not add keyword 'internal' if the
#   class starts with a lower case, e.g. 'matrix'.
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2005-02-11
# o Made all regular expression patterns strict.
# 2005-02-10
# o Rd_parse() in Rdoc$compile() is only called if 'check' is TRUE.
# o Making use of tryCatch() only.
# 2005-02-09
# o If an error occur in tools:Rd_parse() it is now caught and
#   a more informative error is given (with filename and Rdoc name).
# o By default, now Rdoc$compile() runs Rdoc$check() at the end.
# 2005-02-03
# o Added a first simple test for undefined top-level tags in the
#   generated Rd code. Utilizes tools::Rd_parse(), which might be
#   renamed etc. according to its help page.
# 2005-02-02
# o Added support for options() containing tag variables such as @author.
# 2004-10-21
# o Added a local sourceTo() in compile(); will be used in the future.
# 2004-10-18
# o Renamed private argsString() to argsToString().
# o BUG FIX: Rdoc tags was not be parsed for deprecated methods.
# o Added Rdoc comments to all methods.
# 2004-10-13
# o Added check(), which is a wrapper for tools:::check_Rd_files_in_man().
# o setManPath() now calls as.character() too.
# o BUG FIX: setNameFormat("method.class") gave an error.
# o BUG FIX: getPackageNameOf() now only return names of packages and
#   not just any environment.
# o Changed the default name format to "method.class", which is the
#   standard used by R.
# o Updated tagSeemethod() etc to make use of createName().
# 2004-06-27
# o Added \title{} output for RdocDocumention if missing. Similar to
#   RdocAbout.
# 2004-03-11
# o getRdTitle() is now looking for the title in 1) the Rd file in the
#   Rdoc$manPath and then in 2) the CONTENTS files of loaded packages.
# o Added private static getPackageNameOf().
# o getRdMethods() will now also list methods of the given class that are
#   defined in *other* (loaded) packages. This is useful when for
#   instance several packages of the same bundle provides methods to the
#   same class, but one wants to create a list of *all* methods together.
#   This will make it possible to have one bundle called say 'aroma' with
#   packages 'aroma.io', 'aroma.normalize', 'aroma.visualize' and still
#   list all classes in say package 'aroma'.
# 2004-03-10
# o WORK AROUND: Due to a special behavior in LaTeX a verbatim environment
#   may not exist just before \end{Section} without having a line break
#   in between. The problem is that R CMD check & Co removes trailing
#   line breaks. The work around is to add a "\emph{}\n" line to fool it.
#   For Rdoc this is currently done for @howtocite tags. Another reason
#   for using Rdoc ;)
# 2004-03-03
# o Update class to make use of tryCatch() instead of trycatch().
#   This means that R v1.8.1 is now required for R.oo to work correctly.
# 2004-03-02
# o BUG FIX: If variable 'author' did not exists and the @author tag is
#   requested an error was generated, but not the expected RdocException.
#   Instead a missing variable was the problem. Fixed.
# 2003-12-31
# o BUG FIX: For some Rdoc types the \keyword{} statement was placed on
#   the same line as the previous Rd statement. This sometimes generated
#   cluttered Rd index files.
# 2003-11-26
# o Added the tag @RdocDocumentation for general documentation.
# 2003-10-19
# o Now Rdoc tries to create the manPath directory if missing.
# o Added createManPath().
# 2003-09-22
# o Added trial version of @seeothermethods.
# 2003-09-19
# o Now "Methods inherited from Object:" is also sensitive to the
#   showDeprecated argument. By default deprecated methods are not shown.
# o Renamed all arguments 'deprecated' to 'showDeprecated'.
# o Moved createMethods() to Rdoc.DEPRECATED.R.
# 2003-07-18
# o Rdoc$compile() generated an InternalException when a class
#   was not found saying "Not a class". Now it throws an RdocException
#   and is more specific saying that the class does not exists.
#   Updated the Rdoc comments saying pointing out that the classes and
#   methods have to be loaded before calling Rdoc$compile().
# 2003-07-07
# o "BUG FIX": Both @examples and @include does not escape special Rd
#   characters such as '%' -> '\%'.
# o BUG FIX: Forgot to escape '%' in \usage{}, which lead to unbalanced
#   curly brackets when R CMD check * ran.
# 2003-05-03
# o BUG FIX: Misunderstood the usage or \synopsis and \usage. I thought
#   it was either or, but of course not.
# 2003-04-29
# o Rdoc FIX: In the class hierarchy abstract classes was \emph{}'s, but
#   that does not work if they are linked.
# o Now deprecated objects are not generated by default.
# 2003-04-28
# o Added getKeywords() and isKeyword(). @keyword is now checking for
#   existance of keyword.
# o Better error message: The name of source where an error occured is
#   also returned.
# 2003-04-27
# o BUG FIX: Tag values that ended in EOF were not recognized.
# o BUG FIX: tagSet(), tagVisibility() was mistakenly trying to update Rd.
# o BUG FIX: @seemethod did only link to the class, not the method.
# 2003-04-26
# o Major internal update: Much cleaner code and all tags do now rely on
#   function named tag<TagName>().
# 2003-04-26
# o @RdocClass now also generates \keyword{classes}, because promptClass()
#   does it.
# o Update getUsage() to support static methods by returning a string
#   of format "Class$method()".
# 2003-04-23
# o Added get- and setNameFormat(). Now names of the help files can
#   be either "method.class" or "class.method".
# o Improved the Rdoc comments and added a running example.
# 2003-04-13
# o Tags are now ended with a space or a punctuation ( ,;!) but not a
#   period since it can be used in method names and not a colon since
#   that can be used as for instance @see base::help.
# 2003-03-12
# o Added getManPath() and setManPath().
# o Major update: Now Rdoc$compile() will do most of the work. Now the
#   Rdoc perl script basically can be removed and all Rdoc to Rd
#   compilation can be done from within R.
# o Now an <class>.usage.Rdoc file is created for each class.
# o Added private argsString() and public getUsage().
# o Started to write some Rdoc comments for, yes, the Rdoc class.
# 2002-12-08
# o BUG FIX: The generated list of methods in each of the superclasses
#   showed way too many methods. Error was found in private method
#   methodsInheritedFrom().
# 2002-10-21
# o Made Rdoc a class under the new R.oo package.
# 2002-03-05
# * BUG FIX: Added the argument 'this.class=class' in the new() calls in
# both extends() and implements.
# 2002-02-27
# * Found a bug in the code for creating <class>.METHODS. It does not
#   check if the "method" actually is a function. Bug was found in
#   scanForMethods().
# * buildClass() now also creates methods section (<class>.method.Rdoc).
# * Added getTitle() and createMethods().
# 2002-02-26
# * BUG FIX: extends() and implements() did neither work because of
#   interactions with the package "methods"
# * BUG FIX: Due to the new package "methods", my getExtends() did not
#   work as expected in this class. Now getExtends() is only called if
#   the object is of class Object.
# * Updated to make use of setMethodS3().
# 2001-12-29
# * Symbol names for object does not anymore contain a "$" but a ".".
# * Remove require(R.oo). It'll generate a loop in R v 1.4.0.
# 2001-08-06
# * Minor bug fix: Rdoc couldn't generate class documentation for classes
#   in a package that also existed in .GlobalEnv.
# * Made Rdoc a true static class.
# 2001-07-17
# 2001-07-05
# * buildPackage() now loads the library if it is not loaded.
# 2001-07-01
# * Made hierarchy() emphasizes abstract classes, which is more informative.
# * When the perl script Rdoc.pl tried to compile this it complained, since
#   the start and stop tags for Rdoc comments were detected. Made a
#   workaround that works for now.
# 2001-06-29
# * Created!
#########################################################################
