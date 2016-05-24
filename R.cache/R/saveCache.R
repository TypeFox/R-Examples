#########################################################################/**
# @RdocDefault saveCache
#
# @title "Saves data to file cache"
#
# \description{
#  @get "title", which is unique for an optional key object.
# }
#
# @synopsis
#
# \arguments{
#   \item{object}{The object to be saved to file.}
#   \item{key}{An optional object from which a hexadecimal hash
#     code will be generated and appended to the filename.}
#   \item{sources}{Source objects used for comparison of timestamps when
#     cache is loaded later.}
#   \item{suffix}{A @character string to be appended to the end of the
#     filename.}
#   \item{comment}{An optional @character string written in ASCII at the
#     beginning of the file.}
#   \item{pathname}{(Advanced) An optional @character string specifying
#     the pathname to the cache file.  If not specified (default), a unique
#     one is automatically generated from arguments \code{key} and
#     \code{suffix} among other things.}
#   \item{dirs}{A @character @vector constituting the path to the
#      cache subdirectory (of the \emph{cache root directory}
#      as returned by @see "getCacheRootPath") to be used.
#      If @NULL, the path will be the cache root path.}
#   \item{compress}{If @TRUE, the cache file will be saved using
#      gzip compression, otherwise not.}
#   \item{...}{Additional argument passed to @see "base::save".}
# }
#
# \value{
#   Returns (invisible) the pathname of the cache file.
# }
#
# \section{Compression}{
#  The \code{saveCache()} method saves a compressed cache file
#  (with filename extension *.gz) if argument \code{compress} is @TRUE.
#  The @see "loadCache" method locates (via @see "findCache") and
#  loads such cache files as well.
# }
#
# @author
#
# \examples{\dontrun{For an example, see ?loadCache}}
#
# \seealso{
#  For more details on how the hash code is generated etc, @see "loadCache".
# }
#
# @keyword "programming"
# @keyword "IO"
#*/#########################################################################
setMethodS3("saveCache", "default", function(object, key=NULL, sources=NULL, suffix=".Rcache", comment=NULL, pathname=NULL, dirs=NULL, compress=getOption("R.cache::compress", FALSE), ...) {
  # Argument 'compress':
  if (!isTRUE(compress)) compress <- FALSE


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cache file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(pathname)) {
    # Generate cache name from basename and hash object.
    pathname <- generateCache(key=key, suffix=suffix, dirs=dirs);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to file connection
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (compress) {
    pathname <- sprintf("%s.gz", pathname);
    fh <- gzfile(pathname, open="wb");
  } else {
    fh <- file(pathname, open="wb");
  }
  on.exit(close(fh));

  # Save 'identifier'
  identifier <- "Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)";
  if (nchar(identifier) > 64L)
    throw("Internal error. Identifier is too long: ", identifier);
  tail <- paste(rep(" ", times=64L-nchar(identifier)), collapse="");
  identifier <- paste(identifier, tail, sep="");
  writeChar(con=fh, identifier, nchars=64L);

  # Save 'comment'
  if (is.null(comment))
    comment <- "";
  writeBin(con=fh, nchar(comment), size=4L);
  writeChar(comment, con=fh, nchars=nchar(comment));

  # Save 'sources'

  # Look up base::save() once; '::' is expensive
  base_save <- base::save;

  # If 'sources' is not evaluated, it is a so called promise, which will
  # make all of its calling environments to be save too.
  dummy <- is.null(sources);
  base_save(file=fh, sources, compress=compress, ...);

  # Save 'timestamp'
  timestamp <- Sys.time();
  base_save(file=fh, timestamp, compress=compress, ...);

  # Save 'object'
  base_save(file=fh, object, compress=compress, ...);

  invisible(pathname);
})


############################################################################
# HISTORY:
# 2013-12-21
# o Added argument 'pathname' to saveCache().
# 2011-08-16
# o Added support for gzip compressed cache files.
# 2007-01-24
# o Now saveCache() returns the pathname to the cache file.
# 2006-05-25
# o BUG FIX: Work around for not saving "promises" (non-evaluated arguments)
#   in base::save(), which otherwise includes all of the surrounding
#   environment if 'sources' is not evaluated/missing.  For more details
#   see code and my email to r-devel on 2006-05-25.  Thanks to Brian Ripley
#   for explaining what was going on.
# 2006-04-04
# o Added header comment.
# 2005-12-09
# o Object save to file is now a structure containing the object to be
#   cached, a timestamp specifying the it was cached, and a source object.
# o Replaced argument 'file' with 'source'.
# 2005-12-06
# o Created.
############################################################################
