###########################################################################/**
# @RdocFunction .readText
#
# @title "Reads the content of a local or an online text file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @character string specifying a local file or a URL,
#      or a @connection.}
#   \item{...}{Not used.}
#   \item{maxAge}{A @numeric scalar specifying the maximum age (in seconds)
#      of cached URL contents before downloading and recaching.
#      If zero or less, the URL will always be downloaded.}
# }
#
# \value{
#   Returns a @character string.
# }
#
# \section{Caching URL}{
#   When reading online URLs, it may take time a significant time to
#   read its content.  If the content is changing rarely, it is possible
#   to cache the content locally.  If a cached version is found, then it
#   is read instead.
#   It is possible to control how often a file should be recached.  If the
#   cache is older than argument \code{maxAge} (in seconds), then the file
#   is redownloaded and recached.
# }
#
# \section{Newline substitution}{
#   All occurances of \code{\\r\\n} and \code{\\r} are replaced with
#   \code{\\n} such that all lines are ending in \code{\\n} regardless
#   of encoding.
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
.readText <- function(con, ..., maxAge=getOption("R.rsp::downloadIfOlderThan", -Inf)) {
  if (is.character(con)) {
    file <- con;

    # Is the file local and an URL?
    isUrl <- isUrl(file);


    # (a) If URL, download to temporary directory
    if (isUrl) {
      url <- file;
      path <- tempdir();
      filename <- getChecksum(url);
      pathname <- file.path(path, filename);

      # By default, download URL
      download <- TRUE;

      # Unless...
      if (isFile(pathname)) {
        # Age (in seconds) when downloaded file is considered too old
        maxAge <- as.double(maxAge);
        if (is.na(maxAge)) maxAge <- -Inf;
        maxAge <- Arguments$getDouble(maxAge);
        # Time when file was downloaded
        mtime <- file.info(pathname)$mtime;
        # Age of downloaded file in seconds
        dtime <- Sys.time() - mtime;
        units(dtime) <- "secs";
        download <- isTRUE(dtime > maxAge);
      }

      if (download) {
        withoutGString({
          pathname <- downloadFile(url, filename=pathname, skip=FALSE);
        })
      }

      if (isFile(pathname)) file <- pathname;
    } # if (isUrl)


    # (b) Try to open file connection
    con <- tryCatch({
      suppressWarnings({
        file(file, open="rb");
      });
    }, error = function(ex) {
      # (b) If failed, try to download file first
      if (regexpr("^https://", file, ignore.case=TRUE) == -1L) {
        throw(ex);
      }
      url <- file;
      withoutGString({
        pathname <- downloadFile(url, path=tempdir());
      })
      file(pathname, open="rb");
    });
    on.exit(close(con));
  }

  # Sanity check
  stopifnot(inherits(con, "connection"));


  bfr <- NULL;
  while (TRUE) {
    bfrT <- readChar(con, nchars=1e6);
    if (length(bfrT) == 0L) break;
    bfrT <- gsub("\r\n", "\n", bfrT, fixed=TRUE);
    bfrT <- gsub("\r", "\n", bfrT, fixed=TRUE);
    bfr <- c(bfr, bfrT);
  }
  bfr <- paste(bfr, collapse="");
  if (FALSE) {
    bfr <- strsplit(bfr, split="\n", fixed=TRUE);
    bfr <- unlist(bfr, use.names=FALSE);
  }
  bfr;
} # .readText()


##############################################################################
# HISTORY:
# 2013-12-14
# o Added support for caching to .readText().
# o Added Rdoc help.
##############################################################################
