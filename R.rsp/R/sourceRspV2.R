setMethodS3("sourceRspV2", "default", function(..., response=FileRspResponse(file=stdout()), request=NULL, envir=parent.frame(), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'response':
  if (inherits(response, "connection")) {
    response <- FileRspResponse(file=response);
  } else if (is.character(response)) {
    pathname <- Arguments$getWritablePathname(response);
    response <- FileRspResponse(file=pathname);
  } else if (!inherits(response, "RspResponse")) {
    throw("Argument 'response' is not an RspResponse object: ",
                                                         class(response)[1]);
  }

  # Argument 'request':
  if (is.null(request)) {
  } else if (!inherits(request, "HttpRequest")) {
    throw("Argument 'request' is not an HttpRequest object: ", class(request));

  }

  # Argument 'environment':
  if (!is.environment(envir)) {
    throw("Argument 'envir' must be an environment: ", class(envir)[1]);
  }

  # Argument 'verbose':
  if (inherits(verbose, "Verbose")) {
  } else if (is.numeric(verbose)) {
    verbose <- Verbose(threshold=verbose);
  } else {
    verbose <- as.logical(verbose);
    if (verbose)
      verbose <- Verbose(threshold=-1);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # MAIN
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Translate RSP document to an R servlet
  pathnameT <- translateRsp(..., verbose=verbose);
  rCode <- readLines(pathnameT, warn=FALSE);
  attr(rCode, "pathname") <- pathnameT;

  verbose && capture(verbose, paste("22", displayCode(code=rCode)));

  # Parse the servlet from file so that error messages contains line numbers.
  con <- textConnection(rCode);
  on.exit(close(con));

  tryCatch({
    # Parse translated R code
    rExpr <- parse(con);
  }, error = function(ex) {
    # If an parse error occurs, show tranlated code.
    msg <- ex$message;
    line <- gsub(".*line *([0-9]+).*", "\\1", msg);
    code <- displayCode(code=rCode, highlight=line, pager="none");
    code <- unlist(strsplit(code, split="\n", fixed=TRUE), use.names=FALSE);
    ex$code <- code;
    stop(ex);
  })

  # If no 'request' is specifies, create one
  if (is.null(request)) {
    pathname <- attr(rCode, "pathname");
    if (identical(pathname, "")) {
      uri <- "file://stdin/";
    } else if (isUrl(pathname)) {
      uri <- pathname;
    } else {
      uri <- paste("file://", getAbsolutePath(pathname), sep="");
    }
    request <- HttpRequest(uri);
  }

  # Assign 'response' and 'request' to the servlets calling environment.
  assign("response", response, envir=envir);
  assign("request", request, envir=envir);

  # Run servlet
  eval(rExpr, envir=envir);
}, private=TRUE)


##############################################################################
# HISTORY:
# 2011-04-13
# o Added sourceRspV2()
# o Created.
##############################################################################
