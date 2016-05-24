#########################################################################/**
# @RdocDefault readCacheHeader
#
# @title "Loads data from file cache"
#
# \description{
#  @get "title", which is unique for an optional key object.
# }
#
# @synopsis
#
# \arguments{
#   \item{file}{A filename or a @connection.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list structure with element \code{identifier},
#   \code{version}, \code{comment} (optional), \code{sources} (optional),
#   and \code{timestamp}.
# }
#
# @author
#
# @examples "../incl/readCacheHeader.Rex"
#
# \seealso{
#  @see "findCache".
#  @see "loadCache".
#  @see "saveCache".
# }
#
# @keyword "programming"
# @keyword "IO"
# @keyword "internal"
#*/#########################################################################
setMethodS3("readCacheHeader", "default", function(file, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'file':
  if (inherits(file, "connection")) {
    pathname <- sprintf("'%s' [description of the opened connection]", summary(file)$description);
  } else {
    pathname <- as.character(file);
    if (!isFile(pathname))
      throw("Argument 'file' is not an existing file: ", pathname);

    file <- gzfile(pathname, open="rb");
    on.exit({
      if (!is.null(file))
        close(file);
    });
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Try to load cached object from file connection
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  header <- list();

  # 1a. Load identifier:
  id <- readChar(con=file, nchars=64L);
  if (length(id) == 0L) {
    throw("Rcache file format error. Read empty header identifier: ", pathname);
  }
  pattern <- "^Rcache v([0-9][0-9]*[.][0-9][0-9]*([.][0-9][0-9]*)*).*";
  if (regexpr(pattern, id) == -1L) {
    throw("Rcache file format error ('", pathname, "'). Invalid identifier: ", id);
  }
  header$identifier <- id;

  # 1b. Get version
  header$version <- gsub(pattern, "\\1", id);

  # 1c. Read trailing '\0'.
  dummy <- readBin(con=file, what=integer(), size=1, n=1);

  # 2a. Load comment
  if (compareVersion(header$version, "0.1") > 0L) {
    nchars <- readBin(con=file, what=integer(), size=4L, n=1L);
    header$comment <- readChar(con=file, nchars=nchars);
  }

  # 2b. Read trailing '\0'.
  dummy <- readBin(con=file, what=integer(), size=1L, n=1L);

  # 3. Load sources:
  sources <- NULL;  # To please 'codetools' in R v2.6.0
  res <- .baseLoad(con=file, ...);
  header$sources <- res$sources;

  # 4. Load timestamp:
  res <- .baseLoad(con=file, ...);
  header$timestamp <- res$timestamp;

  header;
})


############################################################################
# HISTORY:
# 2012-09-10
# o Updated readCacheHeader() to utilize updated .baseLoad().
# 2011-08-16
# o Added support for reading gzip compressed cache files.
# 2009-10-16
# o Now calling an internal .baseLoad() function of the package.
# 2009-08-11
# o Now readCacheHeader() reports the "pathname" in error/warnings messages,
#   if argument 'file' refers to a file and the "description" if it refers
#   to a connection.
# 2009-07-29
# o ROBUSTNESS: Added sanity check to readCacheHeader() testing that the
#   read header identifier is non-empty.
# 2008-02-14
# o BUG FIX: The throw() for invalid identifiers was trying to put the
#   connection object is the out and not the identifier.
# 2006-04-04
# o Created.
############################################################################
