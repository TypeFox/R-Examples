###########################################################################/**
# @RdocDefault browseRsp
# @alias browseRsp.Package
#
# @title "Starts the internal web browser and opens the URL in the default web browser"
#
# \description{
#  @get "title".
#  From this page you access not only help pages and demos on how to use
#  RSP, but also other package RSP pages.
# }
#
# @synopsis
#
# \arguments{
#   \item{url}{A @character string for the URL to be viewed.
#     By default the URL is constructed from the \code{urlRoot} and
#     the \code{path} parameters.
#   }
#   \item{urlRoot}{A @character string specifying the URL root.  By default
#     the URL is constructed from the \code{host} and the \code{port}.}
#   \item{host}{An optional @character string for the host of the URL.}
#   \item{port}{An optional @integer for the port of the URL.}
#   \item{path}{An optional @character string for the context path of the URL.}
#   \item{start}{If @TRUE, the internal \R web server is started if not
#     already started, otherwise not.}
#   \item{stop}{If @TRUE, the internal \R web server is stopped, if started.}
#   \item{...}{Additional arguments passed to @see "utils::browseURL".}
# }
#
# \value{
#   Returns (invisibly) the URL.
# }
#
# @author
#
# \seealso{
#   Internally, @see "utils::browseURL" is used to launch the browser.
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("browseRsp", "default", function(url=paste(urlRoot, path, sep="/"), urlRoot=sprintf("http://%s:%d", host, port), host="127.0.0.1", port=8074L, path="", start=TRUE, stop=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'url':
  if (!is.null(url)) {
    if (!isUrl(url) && isUrl(urlRoot) && nchar(path) == 0L) {
      path <- getRelativePath(url);
      if (isAbsolutePath(path)) {
        throw("Cannot open file, because it is not possible to infer its relative pathname: ", url);
      }
      url <- paste(urlRoot, path, sep="/");
    }
  }

  # Argument 'port':
  port <- Arguments$getInteger(port, range=c(0L,65535L));

  # Argument 'stop':
  stop <- Arguments$getLogical(stop);


  # Get the/a HTTP daemon
  httpDaemon <- getStaticInstance(HttpDaemon);

  # Stop HTTP server?
  if (stop) {
    if (isStarted(httpDaemon))
      stop(httpDaemon);
    return(!isStarted(httpDaemon));
  }

  # Start HTTP server?
  if (start) {
    # Add the following directories to the list of known root paths:
    # (1) current directory
    paths <- getwd();
    # (2) rsp/ under current directory
    paths <- c(paths, file.path(getwd(), "rsp"));
    # (3) the parent of all library paths
    paths <- c(paths, dirname(.libPaths()));
    # (4) /library/R.rsp/rsp/ as a root path too.
    paths <- c(paths, system.file("rsp", package="R.rsp"));
    appendRootPaths(httpDaemon, paths);

    if (!isStarted(httpDaemon)) {
      # Start the web server
      start(httpDaemon, port=port, default="^index[.](html|.*)$")
    }
  }

  if (!is.null(url)) {
    browseURL(url, ...);
  }

  invisible(url);
})


setMethodS3("browseRsp", "Package", function(this, ..., path=sprintf("library/%s/rsp/", getName(this))) {
  browseRsp(..., path=path);
})


############################################################################
# HISTORY:
# 2013-05-22
# o Now it is possible to call R.rsp::browseRsp() without loading R.rsp.
# o Now browseRsp() automatically adds the current directory as well
#   an rsp/ subdirectory of it to the list of known root paths.
# o Now browseRsp("foo.html") is interpreted as browseRsp(path="foo.html")
#   whenever then requested file is not in an URL format.
# o Added argument 'urlRoot' to browseRsp().
# 2011-04-18
# o Now browseRsp() pass '...' arguments to browseURL().
# 2011-03-12
# o Replaced all references to static class HttpDaemon by a single one.
#   This will make it easier to generalize the code in the future.
# 2007-07-19
# o Added Rdoc comments.
# 2007-07-11
# o Added browseRsp() for the Package too, e.g. browseRsp(aroma.light).
# o Now (the parents of) all library paths are added to the root paths.
# o Now the root paths are updated each time browseRsp(start=TRUE) is
#   called, even if the HTTP daemon is already started.
# 2005-10-18
# o Created.
############################################################################
