###########################################################################/**
# @RdocClass HttpDaemon
#
# @title "The HttpDaemon class"
#
# \description{
#  @classhierarchy
#
#  A minimalistic HTTP daemon (web server) that also preprocesses RSP.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# \details{
#  The actual server is written in Tcl such that it runs in a non-blocking
#  mode, which means that the R prompt will be available for other things.
#  This class is tightly coupled with the source code of the Tcl script.
#
#  For security reasons, the server only accept connections from the
#  local host (127.0.0.1).  This lowers the risk for external computers
#  to gain access to the R session.
#  This is asserted by the \code{accept_connect} Tcl procedure in
#  r-httpd.tcl (located in \code{system("tcl/", package="R.rsp")}).
#  If access from other hosts are wanted, then this procedure needs to
#  be modified.
#
#  The Tcl server was written by Steve Uhlers, and later adopted for R by
#  Philippe Grosjean and Tom Short (Rpad package author) [1].
# }
#
# @examples "../incl/HttpDaemon.Rex"
#
# \references{
#   [1] Rpad package, Tom Short, 2005.\cr
# }
#
# @author
#
# @keyword IO
# @keyword internal
#*/###########################################################################
setConstructorS3("HttpDaemon", function(...) {
  this <- extend(Object(), "HttpDaemon",
    .debug = FALSE,
    .count = 0L,
    .rootPaths = NULL
  )

  this$count <- this$count + 1L;

  # Check if another server is already running.
  if (this$count > 1L) {
    throw("ERROR: There is already an HttpDaemon running. Sorry, but the current implementation allows only one per R session.");
  }

  this;
})

setMethodS3("finalize", "HttpDaemon", function(this, ...) {
  if (isStarted(this))
    stop(this);
  this$count <- this$count - 1L;
}, protected=TRUE, createGeneric=FALSE)


setMethodS3("getCount", "HttpDaemon", function(static, ...) {
  as.integer(static$.count);
}, protected=TRUE)


setMethodS3("setCount", "HttpDaemon", function(static, count, ...) {
  static$.count <- as.integer(count);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the HTTP daemon"
#
# \description{
#  @get "title".
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
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "HttpDaemon", function(x, ...) {
  # To please R CMD check
  static <- x;

  s <- paste(class(static)[1L], ":", sep="");
  if (isStarted(static)) {
    s <- paste(s, " HTTP daemon is started.", sep="");
    s <- paste(s, " Current root paths: ", paste(getRootPaths(static), collapse=";"), ".", sep="");
    s <- paste(s, " Port: ", getPort(static), ".", sep="");
    s <- paste(s, " Default filename: ", getDefaultFilenamePattern(static),
                                                        ".", sep="");
  } else {
    s <- paste(s, " HTTP daemon is not started.", sep="");
  }
  s;
})



#########################################################################/**
# @RdocMethod openUrl
#
# @title "Starts the HTTP daemon and launches the specified URL"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{url}{The URL to be opened.}
#   \item{host}{The host where the HTTP server is running.}
#   \item{port}{The port to be used.}
#   \item{path}{The path to the document to be opened.}
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
#   Called by for instance @seemethod "startHelp".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("openUrl", "HttpDaemon", function(static, url=sprintf("http://%s:%d/%s", host, port, path), host="127.0.0.1", port=8074, path="", ...) {
  # - - - - - - - g- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'port':
  port <- Arguments$getInteger(port, range=c(0,65535));


  # Start HTTP server, if not started.
  if (!isStarted(static)) {
    # Start the web server
    rootPath <- system.file("rsp", package="R.rsp")
    start(static, rootPath=rootPath, port=port, ...);
  }

  if (!is.null(url))
    browseURL(url);
})


#########################################################################/**
# @RdocMethod startHelp
#
# @title "Starts the HTTP daemon and launches the help page"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @seemethod "openUrl".}
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
# @keyword IO
#*/#########################################################################
setMethodS3("startHelp", "HttpDaemon", function(static, ...) {
  openUrl(static, path="R/Help/", ...);
})








#########################################################################/**
# @RdocMethod getConfig
#
# @title "Retrieves the server's 'config' structure from Tcl"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a tclArray object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getConfig", "HttpDaemon", function(static, ...) {
  # Load required package
  requireNamespace("tcltk") || stop("Package not installed/found: tcltk");

  config <- tcltk::as.tclObj("config");
  class(config) <- c("tclArray", class(config));
  config;
}, static=TRUE, protected=TRUE)





#########################################################################/**
# @RdocMethod getHttpRequest
#
# @title "Gets the HTTP request"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "HttpRequest" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getHttpRequest", "HttpDaemon", function(static, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getRequestUri <- function(...) {
    url <- NA;
    tryCatch({
      url <- as.character(tcltk::tclvalue("url"));
    }, error = function(ex) {
    })
    url;
  }

  getData <- function(field=NULL, ...) {
    data <- tcltk::as.tclObj("data");
    class(data) <- c("tclArray", class(data));
    if (is.null(field))
      return(data);
    value <- data[[field]];
    if (is.null(value))
      return(NULL);
    value <- tcltk::tclvalue(value);
    value;
  }

  getRequestParameters <- function(...) {
    params <- list();
    query <- getData("query");
    if (!is.null(query)) {
      query <- strsplit(query, split="&", fixed=TRUE)[[1L]];
      if (length(query) == 0L)
        return(params);

      query <- strsplit(query, split="=", fixed=TRUE);

      for (kk in seq(along=query)) {
        pair <- query[[kk]];
        name <- urlDecode(pair[1L]);
        value <- urlDecode(pair[2L]);
        params[[kk]] <- value;
        names(params)[kk] <- name;
      }
    }

    params;
  }

  HttpRequest(
    serverPort    = getPort(static),
    contextRoot   = getParent(as.character(tcltk::tclvalue("mypath"))),
    requestUri    = getData("url"),
    queryString   = getData("query"),
    remoteAddress = getData("ipaddr"),
    parameters    = getRequestParameters(static)
  )
}, static=TRUE)




#########################################################################/**
# @RdocMethod getPort
#
# @title "Gets the socket port of the HTTP daemon"
#
# \description{
#  @get "title", if started.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer if started, otherwise @NA.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getPort", "HttpDaemon", function(static, ...) {
  config <- getConfig(static);
  as.integer(config$port);
}, static=TRUE)




#########################################################################/**
# @RdocMethod getRootPaths
#
# @title "Gets the root directories of the HTTP daemon"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @vector of @character string if started, otherwise @NA.
# }
#
# @author
#
# \seealso{
#   @seemethod setRootPaths
#   @seemethod appendRootPaths
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getRootPaths", "HttpDaemon", function(static, ...) {
  # If server is started, updated rootPaths from the servers settings
  if (isStarted(static)) {
    paths <- tcltk::tcl("getRootPaths");
    paths <- as.character(paths);
    static$.rootPaths <- paths;
  }

  static$.rootPaths;
}, static=TRUE)




#########################################################################/**
# @RdocMethod setRootPaths
#
# @title "Sets a new set of root directories for the HTTP daemon"
#
# \description{
#  @get "title", if started.
# }
#
# @synopsis
#
# \arguments{
#   \item{paths}{A @vector of paths.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns (invisibly) the previously known root directories.
# }
#
# @author
#
# \seealso{
#   @seemethod getRootPaths
#   @seemethod appendRootPaths
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("setRootPaths", "HttpDaemon", function(static, paths, ...) {
  oldPaths <- getRootPaths(static);

  # Keep only unique paths
  paths <- unlist(strsplit(paths, split=";", fixed=TRUE), use.names=FALSE);
  paths <- unique(paths);
  static$.rootPaths <- paths;

  # If server is started, updated servers settings
  if (isStarted(static)) {
    paths <- paste(paths, collapse=";");
    res <- tcltk::tcl("setRootPaths", paths);
  }

  invisible(oldPaths);
}, static=TRUE)


## setMethodS3("refreshRootPaths", "HttpDaemon", function(static, ...) {
##   # If server is started, updated servers settings
##   if (isStarted(static)) {
##     paths <- getRootPaths(static);
##     paths <- paste(paths, collapse=";");
##     res <- tcltk::tcl("setRootPaths", paths);
##   }
##   invisible(getRootPaths(static));
## }, static=TRUE)



#########################################################################/**
# @RdocMethod appendRootPaths
# @aliasmethod insertRootPaths
#
# @title "Appends and inserts new paths to the list of known root directories"
#
# \description{
#  @get "title", if started.
# }
#
# @synopsis
#
# \arguments{
#   \item{paths}{A @vector of paths.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns (invisibly) the previously known root directories.
# }
#
# @author
#
# \seealso{
#   @seemethod getRootPaths
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("appendRootPaths", "HttpDaemon", function(static, paths, ...) {
  setRootPaths(static, c(getRootPaths(static), paths), ...);
}, static=TRUE)


setMethodS3("insertRootPaths", "HttpDaemon", function(static, paths, ...) {
  setRootPaths(static, c(paths, getRootPaths(static)), ...);
}, static=TRUE)





#########################################################################/**
# @RdocMethod getDefaultFilenamePattern
#
# @title "Gets the default filename pattern to be loaded by the HTTP daemon"
#
# \description{
#  @get "title", if started.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @character string if started, otherwise @NA.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getDefaultFilenamePattern", "HttpDaemon", function(static, ...) {
  config <- getConfig(static);
  as.character(config$default);
}, static=TRUE)




#########################################################################/**
# @RdocMethod isStarted
#
# @title "Checks if the HTTP daemon is started"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns @TRUE if the server is started, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seemethod "start" and @seemethod "stop".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("isStarted", "HttpDaemon", function(x, ...) {
  # To please R CMD check...
  static <- x;

  port <- getPort(static);
  (length(port) != 0L);
}, static=TRUE)




#########################################################################/**
# @RdocMethod sourceTcl
#
# @title "Loads the Tcl source for the HTTP daemon into R"
#
# \description{
#  @get "title".
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
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("sourceTcl", "HttpDaemon", function(static, ...) {
  # Load required package
  requireNamespace("tcltk") || stop("Package not installed/found: tcltk");

  tclPath <- system.file("tcl", package="R.rsp");
  pathname <- file.path(tclPath, "r-httpd.tcl");
  if (!isFile(pathname))
    stop("Tcl source code file not found: ", pathname);

  res <- tcltk::tcl("source", pathname);
  invisible(res);
}, static=TRUE, protected=TRUE);




#########################################################################/**
# @RdocMethod start
#
# @title "Starts the HTTP daemon"
#
# \description{
#  @get "title".  Currently, only one HTTP daemon can run at each time,
#  regardless of port used.
# }
#
# @synopsis
#
# \arguments{
#   \item{rootPaths}{The path(s) to act as the root of the webserver file
#       system.  Files in parent directories of the root, will not be
#       accessable.  If @NULL, the preset paths will be used,
#       cf. @seemethod "setRootPaths".}
#   \item{port}{The socket port the server listens to.}
#   \item{default}{The default filename pattern to be retrieved if
#       not specified.}
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
#   @seemethod "setRootPaths".
#   @seemethod "isStarted".
#   @seemethod "stop".
#   @seemethod "restart".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("start", "HttpDaemon", function(x, rootPaths=NULL, port=8080, default="^index[.](html|.*)$", ...) {
  # The R.rsp package needs to be attached in order to make certain
  # R functions of R.rsp available to the Tcl HTTP daemon.
  use("R.rsp", quietly=TRUE)

  # To please R CMD check...
  static <- x;

  # Is HTTP daemon already started?
  if (isStarted(static))
    stop("HTTP daemon is already started: ", as.character(static));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rootPaths':
  if (length(rootPaths) > 0L) {
    rootPaths <- unlist(strsplit(rootPaths, split=";", fixed=TRUE), use.names=FALSE);
    rootPaths <- unlist(sapply(rootPaths, FUN=function(path) {
      Arguments$getReadablePathname(path, mustExist=TRUE);
    }), use.names=FALSE)
    setRootPaths(static, rootPaths);
  } else {
    rootPaths <- getRootPaths(static);
  }

  # Argument 'port':
  port <- Arguments$getInteger(port, range=c(0L,65535L));

  # Argument 'default':
  default <- Arguments$getCharacter(default, nchar=c(1L,256L));

  # Source the TCL httpd code
  sourceTcl(static);

  # Start the HTTP daemon (the webserver)
  res <- tcltk::tcl("server", paste(rootPaths, collapse=";"), port, default);

  # Validate opened port.
  port <- Arguments$getInteger(tcltk::tclvalue(res), range=c(0L,65535L));

  invisible(port);
}, static=TRUE, createGeneric=FALSE)




#########################################################################/**
# @RdocMethod stop
#
# @title "Stops the HTTP daemon"
#
# \description{
#  @get "title".
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
#   @seemethod "isStarted".
#   @seemethod "start".
#   @seemethod "restart".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("stop", "HttpDaemon", function(static, ...) {
  # Is HTTP daemon already started?
  if (!isStarted(static))
    stop("HTTP daemon is not started.");

  # Close the httpd socket.
  tcltk::.Tcl("close $config(listen)");
  tcltk::.Tcl("unset config");

  invisible(TRUE);
}, static=TRUE)




#########################################################################/**
# @RdocMethod restart
#
# @title "Restarts the HTTP daemon"
#
# \description{
#  @get "title".
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
#   @seemethod "isStarted".
#   @seemethod "start".
#   @seemethod "stop".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("restart", "HttpDaemon", function(static, ...) {
  if (!isStarted(static))
    throw("HTTP daemon not started.");

  rootPaths <- getRootPaths(static);
  port <- getPort(static);
  default <- getDefaultFilenamePattern(static);

  stop(static, ...);

  start(static, rootPaths=rootPaths, port=port, default=default, ...);
}, static=TRUE)




#########################################################################/**
# @RdocMethod writeResponse
#
# @title "Writes a string to the HTTP output connection"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{A set of @character strings to be outputted.}
# }
#
# \details{
#   \emph{Note: For efficiency, there is no check if the HTTP daemon is
#         started or not.}
# }
#
# \value{
#  Returns (invisibly) the number of characters written.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("writeResponse", "HttpDaemon", function(static, ...) {
  str <- paste(..., collapse="", sep="");

  # Nothing to do?
  if (nchar(str) == 0L) {
    return(invisible(0L));
  }

  if (isTRUE(static$.debug)) {
    mcat("=========================================================\n");
    mcat("= BEGIN: Fake HttpDaemon response\n");
    mcat("=========================================================\n");
    mcat(str);
    mcat("=========================================================\n");
    mcat("= END: Fake HttpDaemon response\n");
    mcat("=========================================================\n");
  } else {
    # Escape certain characters, by converting the string to a Tcl string
    # and back.
    str <- as.character(tcltk::tclVar(str));

    # Write the string to HTTP output connection.
    tcltk::.Tcl(paste("catch { puts $sock $", str, " }", sep=""));
  }

  invisible(nchar(str));
})


###############################################################################
# HISTORY:
# 2013-09-18
# o ROBUSTNESS: Now start() for HttpDaemon makes sure that the R.rsp package
#   is attached so that the Tcl HTTP daemon have access to its methods.
# 2013-03-31
# o Now HttpDaemon$openUrl() passes '...' to start().
# o Now HttpDaemon$start() uses default="^index[.](html|.*)$".
# o Renamed getDefaultFilename() to getDefaultFilenamePattern().
# 2013-02-23
# o Now writeResponse() for HttpDaemon writes to standard output only,
#   if HttpDaemon$.fake is TRUE.
# 2011-09-21
# o BUG FIX: HttpDaemon$getRootPaths() did not handle paths with
#   spaces correctly.  Added a getRootPaths() Tcl function to
#   instead handle this, which is called by the former.
# 2011-03-12
# o CLEANUP: Replaced on HttpDaemon$<method>(...) with <method>(static, ...).
# 2011-03-08
# o BUG FIX: getHttpRequest() for HttpDaemon would drop all but the last
#   of replicated query parameters of the same name.  Thanks to Truc Trung
#   at University of Bergen, Norway, for reporting on this.
# 2011-01-06
# o DOCUMENTATION: Clarified in the help of HttpDaemon that it is only
#   connections from the local host (127.0.0.1) that are accepted.
#   This lowers the risk for unauthorized access to the R session.
# 2007-06-10
# o Now all methods of 'tcltk' are called explicitly with prefix 'tcltk::'.
# 2006-07-10
# o Now append- and insertRootPaths() pass arguments '...' to setRootPaths().
# 2006-07-04
# o Added openUrl().
# 2006-10-13
# o BUG FIX: Used obsolete setClassS3() instead of setConstructorS3().
# 2006-01-21
# o Added writeResponse().
# o Moved processRsp() to its own file.  The purpose is to one day get a
#   HttpDaemon class which does not know of RSP pages.
# 2006-01-12
# o The example of HttpDaemon is now "runnable" in interactive mode.
# 2005-11-30
# o Added restart().
# o Now processRsp() uses new HttpDaemonResponse class which outputs written
#   response directly to the Tcl HTTP Daemon output stream.  This is one step
#   closer to a immediate output to the browser.
# 2005-10-20
# o Now root paths can be set before the server has started.
# 2005-10-19
# o Added append- and insertRootPaths().  Modified the server so it supports
#   multiple root directories.
# o Now the search for an existing file is 100% done by the Tcl HTTP daemon.
# 2005-09-26
# o Added Rdoc comments to getHttpRequest().
# o Removed getRequestContext() and getRequestParameters().
# 2005-09-24
# o Added getRequestContext(), getRequestParameters(), getRequest().
# 2005-09-22
# o Added Rdoc comments.
# o Added RSP preprocessor. It really works! Sweet.
# o Created static HttpDaemon class.
# o Added rsp to list of known mimetypes. /HB
# o Adopted from the minihttpd.tcl in the Rpad package. /HB
###############################################################################
