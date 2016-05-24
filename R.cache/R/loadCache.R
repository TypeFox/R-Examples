#########################################################################/**
# @RdocDefault loadCache
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
#   \item{key}{An optional object from which a hexadecimal hash
#     code will be generated and appended to the filename.}
#   \item{sources}{Optional source objects.  If the cache object has a
#     timestamp older than one of the source objects, it will be ignored
#     and removed.}
#   \item{suffix}{A @character string to be appended to the end of the
#     filename.}
#   \item{removeOldCache}{If @TRUE and the cache is older than the
#     \code{sources}, the cache file is removed, otherwise not.}
#   \item{pathname}{The pathname to the cache file.  If specified,
#     arguments \code{key} and \code{suffix} are ignored.  Note that
#     this is only needed in order to read a cache file for which
#     the key is unknown, for instance, in order to investigate
#     an unknown cache file.}
#   \item{dirs}{A @character @vector constituting the path to the
#      cache subdirectory (of the \emph{cache root directory}
#      as returned by @see "getCacheRootPath") to be used.
#      If @NULL, the path will be the cache root path.}
#   \item{...}{Not used.}
#   \item{onError}{A @character string specifying what the action is
#      if an exception is thrown.}
# }
#
# \value{
#   Returns an \R object or @NULL, if cache does not exist.
# }
#
# \details{
#   The hash code calculated from the \code{key} object is a
#   32 characters long hexadecimal MD5 hash code.
#   For more details, see @see "getChecksum".
# }
#
# @author
#
# @examples "../incl/loadCache.Rex"
#
# \seealso{
#  @see "saveCache".
# }
#
# @keyword "programming"
# @keyword "IO"
#*/#########################################################################
setMethodS3("loadCache", "default", function(key=NULL, sources=NULL, suffix=".Rcache", removeOldCache=TRUE, pathname=NULL, dirs=NULL, ..., onError=c("warning", "print", "quiet", "error")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'onError':
  onError <- match.arg(onError);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find cached file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(pathname))
    pathname <- findCache(key=key, suffix=suffix, dirs=dirs);
  if (is.null(pathname))
    return(NULL);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Try to load cached object from file connection
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!isFile(pathname))
    return(NULL);

  fh <- gzfile(pathname, open="rb");
  on.exit({
    if (!is.null(fh))
      close(fh);
  });

  tryCatch({
    header <- readCacheHeader(fh);
    if (!is.null(sources))
      header$sources <- sources;

    timestamp <- NULL;
    attachLocally(header);  # Attaches 'timestamp' (and 'sources')

    if (!is.null(timestamp) && !is.null(sources)) {
      for (sourcePathname in sources) {
        if (!is.character(sources)) {
          warning("No timestamp check of cache was performed. Unsupported type of cache source: ", class(sources)[1]);
          break;
        }

        if (!file.exists(sourcePathname)) {
          warning("No timestamp check of cache was performed. Source file not found: ", sourcePathname);
          break;
        }

        # Is source file newer than cache?
        lastModified <- file.info(sourcePathname)$mtime;
        if (lastModified > timestamp) {
          # Remove out-of-date cache file?
          if (removeOldCache) {
            close(fh); fh <- NULL;
            file.remove(pathname);
          }
          return(NULL);
        }
      } # for (sourcePathname in sources)
    }

    # 4. Load cached object:
    res <- .baseLoad(con=fh, ...);
    object <- res$object;
    res <- NULL; # Not needed anymore

    # 5. Update the "last-modified" timestamp of the cache file?
    touch <- getOption("R.cache::touchOnLoad");
    touch <- identical(touch, TRUE);
    if (touch) {
      touchFile(pathname);
    }

    # 6. Return cached object
    return(object);
  }, error = function(ex) {
     if (onError == "print") {
       print(ex);
     } else if (onError == "warning") {
       warning(ex);
     } else if (onError == "error") {
       stop(ex);
     }
  })

  NULL;
})


############################################################################
# HISTORY:
# 2012-09-10
# o Updated readCacheHeader() to utilize updated .baseLoad().
# 2011-08-16
# o Added support for loading gzip compressed cache files.
# 2009-10-16
# o Now calling an internal .baseLoad() function of the package.
# 2009-09-11
# o Added argument 'onError' to loadCache(), to specify the action when
#   an error occurs.  The default used to be to print the error message
#   (onError="print"), but now the default is to generate a warning
#   ("warning").  The other alternatives are do silently ignore it, or
#   to throw the error ("error").  Except for onError="error", loadCache()
#   always returns NULL if an error occurs.
# 2008-02-27
# o Added option to updated the "last-modified" timestamp of cache files
#   whenever they are loaded via loadCache().
# 2008-02-14
# o Now errors reports the pathname, if available.
# 2006-08-09
# o Added link to cache() in Biobase.
# 2006-05-25
# o Added argument 'pathname' to make it possible to load a cache file
#   explicitly.
# 2006-04-04
# o Added header comment for file format > v0.1.
# o Added detection of file format version.
# 2005-12-09
# o Added support for internal 'cache' and 'timestamp' objects.
# 2005-12-06
# o Created.
############################################################################
