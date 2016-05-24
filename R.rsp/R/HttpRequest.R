###########################################################################/**
# @RdocClass HttpRequest
#
# @title "The HttpRequest class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{requestUri}{A @character string of the requested URI.}
#   \item{parameters}{A named @list of parameter values.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
# @keyword internal
#*/###########################################################################
setConstructorS3("HttpRequest", function(requestUri=NULL, parameters=list(), ...) {
  if (is.list(requestUri)) {
    request <- requestUri;
    requestUri <- request$requestUri;
    parameters <- request$parameters;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'requestUri':
  requestUri <- Arguments$getCharacter(requestUri);

  # Argument 'parameters':
  if (!is.list(parameters))
    stop("Argument 'parameters' must be a list: ", mode(parameters));

  extend(Object(), "HttpRequest",
    serverPort = NA,
    serverName = NA,
    contextRoot = ".",
    contextType = NA,
    contextLength = -1,
    remoteAddress = NA,
    remoteHost = NA,
    scheme = NA,
    protocol = NA,
    requestUri = requestUri,
    parameters = parameters,
    ...
  )
})


###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the HTTP request"
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
setMethodS3("as.character", "HttpRequest", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- paste(class(this)[1], ":", sep="");

  if (is.null(this$requestUri)) {
    s <- paste(s, " Request URI: <none>.", sep="");
  } else {
    s <- paste(s, " Request URI: ", this$requestUri, ".", sep="");
  }

  if (nbrOfParameters(this) > 0) {
    params <- unlist(this$parameters, use.names=TRUE);
    params <- paste(names(params), params, sep="=");
    params <- paste(params, collapse=", ");
    s <- paste(s, " Parameters: ", params, ".", sep="");
  } else {
    s <- paste(s, " Parameters: <none>.", sep="");
  }
  s;
})





#########################################################################/**
# @RdocMethod nbrOfParameters
#
# @title "Gets the number of parameters"
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
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seemethod "getParameter".
#   @seemethod "hasParameter".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("nbrOfParameters", "HttpRequest", function(this, ...) {
  length(this$parameters);
})



#########################################################################/**
# @RdocMethod getParameters
#
# @title "Gets all parameters"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{trim}{If @TRUE, each parameter value is trimmed of whitespace.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a named @list.
# }
#
# @author
#
# \seealso{
#   @seemethod "getParameter".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getParameters", "HttpRequest", function(this, trim=FALSE, ...) {
  params <- as.list(this$parameters);
  if (trim) {
    params <- lapply(params, FUN=trim);
  }
  params;
})



#########################################################################/**
# @RdocMethod getParameter
#
# @title "Gets a parameter"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{Name of parameter to be retrieved.}
#   \item{default}{Value to be returned if parameter is missing.}
#   \item{drop}{If @TRUE and the number of returned values is one, then
#    this single value is returned, otherwise a named @vector.}
#   \item{...}{Additional arguments passed to @seemethod "getParameters".}
# }
#
# \value{
#  Returns the value(s) of the parameter either as a single value or
#  as a named @list.
#  If the parameter does not exist, the default value is returned as is.
# }
#
# @author
#
# \seealso{
#   @seemethod "hasParameter".
#   @seemethod "getParameters".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getParameter", "HttpRequest", function(this, name, default=NULL, drop=TRUE, ...) {
  if (hasParameter(this, name)) {
    params <- getParameters(this, ...);
    idxs <- which(names(params) == name);
    params <- params[idxs];

    if (drop && length(params) == 1L) {
      params <- params[[1L]];
    }
  } else {
    params <- default;
  }

  params;
})





#########################################################################/**
# @RdocMethod hasParameter
#
# @title "Checks if a parameter exists"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{Name of parameter to be checked.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns @TRUE if the parameter exists, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seemethod "getParameter".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("hasParameter", "HttpRequest", function(this, name, ...) {
  name <- Arguments$getCharacter(name, nchar=c(1,256));
  is.element(name, names(this$parameters));
})





#########################################################################/**
# @RdocMethod getRemoteAddress
#
# @title "Gets the IP address of the client that sent the request"
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
#   @see "getRemoteHost".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getRemoteAddress", "HttpRequest", function(this, ...) {
  this$remoteAddress;
})



#########################################################################/**
# @RdocMethod getRemoteHost
#
# @title "Gets the fully qualified name of the client that sent the request"
#
# \description{
#  @get "title".
#  If it cannot resolve the hostname, this method returns the dotted-string
#  form of the IP address.
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
#   @see "getRemoteAddress".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getRemoteHost", "HttpRequest", function(this, ...) {
  this$remoteHost;
})




#########################################################################/**
# @RdocMethod getServerName
#
# @title "Gets the host name of the server that revieved the request"
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
#   @see "getServerPort".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getServerName", "HttpRequest", function(this, ...) {
  this$serverName;
})




#########################################################################/**
# @RdocMethod getServerPort
#
# @title "Gets the port number on which this request was received"
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
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @see "getServerPort".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getServerPort", "HttpRequest", function(this, ...) {
  as.integer(this$serverPort);
})



#########################################################################/**
# @RdocMethod getScheme
#
# @title "Gets the scheme used to make this request"
#
# \description{
#  @get "title", e.g. http, https, or ftp.
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
#   @see "getProtocol".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getScheme", "HttpRequest", function(this, ...) {
  this$scheme;
})


#########################################################################/**
# @RdocMethod getProtocol
#
# @title "Gets the name and version of the protocol used to make this request"
#
# \description{
#  @get "title", e.g. HTTP/1.1.
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
#   @see "getScheme".
#   @seeclass
# }
#
# @keyword IO
#*/#########################################################################
setMethodS3("getProtocol", "HttpRequest", function(this, ...) {
  this$protocol;
})



#########################################################################/**
# @RdocMethod getContentType
#
# @title "Gets the MIME type of the body of the request"
#
# \description{
#  @get "title", or @NULL if the type is not known.
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
#*/#########################################################################
setMethodS3("getContentType", "HttpRequest", function(this, ...) {
  this$contentType;
})


#########################################################################/**
# @RdocMethod getContentLength
#
# @title "Gets the length of contents"
#
# \description{
#  @get "title" (in bytes), or -1 if the length is not known.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer.
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
setMethodS3("getContentLength", "HttpRequest", function(this, ...) {
  len <- this$contentLength;
  if (is.null(len))
    len <- -1;
  as.integer(len);
})

setMethodS3("getDateHeader", "HttpRequest", function(this, ...) {
}, protected=TRUE)

setMethodS3("getHeader", "HttpRequest", function(this, ...) {
}, protected=TRUE)


setMethodS3("getContextPath", "HttpRequest", function(this, ...) {
}, protected=TRUE)


setMethodS3("getQueryString", "HttpRequest", function(this, ...) {
  this$queryString;
}, protected=TRUE)

setMethodS3("getRemoteUser", "HttpRequest", function(this, ...) {
}, protected=TRUE)

setMethodS3("getRequestUri", "HttpRequest", function(this, ...) {
  this$requestUri;
}, protected=TRUE)

setMethodS3("getRequestUrl", "HttpRequest", function(this, ...) {
}, protected=TRUE)

setMethodS3("getServletPath", "HttpRequest", function(this, ...) {
}, protected=TRUE)



#########################################################################/**
# @RdocMethod getRealPath
#
# @title "Gets the file system path for a given URI"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{uri}{A URI as a @character string.}
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
#*/#########################################################################
setMethodS3("getRealPath", "HttpRequest", function(this, uri, ...) {
  contextRoot <- this$contextRoot;
  realPath <- filePath(contextRoot, uri);
  realPath;
})


##############################################################################
# HISTORY:
# 2013-05-22
# o Added argument 'trim=FALSE' to getParameter() and getParameters()
#   for HttpRequest.
# 2011-03-08
# o Updated getParameter() of HttpRequest to returning the value of a
#   query parameters with multiple entries.  Added argument 'drop'.
# 2006-02-22
# o Added getParameters() for completeness.
# 2005-10-27
# o Added missing Rdoc comments.
# 2005-09-26
# o Added several getNNN() methods.
# 2005-09-24
# o Created.
#############################################################################
