###########################################################################/**
# @RdocClass Exception
#
# \title{The Exception class to be thrown and caught}
#
# \description{
#  @classhierarchy
#
#  Creates an Exception that can be thrown and caught. The \code{Exception}
#  class is the root class of all other \code{Exception} classes.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{One or several strings, which will be concatenated and contain
#     informative message about the exception.}
#   \item{sep}{The string to used for concatenating several strings.}
#   \item{collapse}{The string to used collapse vectors together.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @examples "../incl/Exception.Rex"
#
# @author
#
# \seealso{
#   See also \code{\link[base:conditions]{tryCatch}()} (and @see "base::try").
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setConstructorS3("Exception", function(..., sep="", collapse=", ") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  sys.functions <- function(parents) {
    functions <- list();
    for (kk in seq_along(parents)) {
      parent <- parents[kk];
      functions[[kk]] <- sys.function(which=kk);
    }
    functions;
  } # sys.functions()

  fcnPathname <- function(call) {
    res <- attr(call, "srcref");
    if (is.null(res)) return("");
    res <- attr(res, "srcfile");
    res$filename;
  } # fcnPathname()

  fcnName <- function(call) {
    code <- deparse(call[1]);
#    code <- grep("^function\\(", code, value=TRUE);
    if (length(code) == 0) return("");
    code <- code[1];
    code <- gsub("\\(.*", "", code);
    code;
  } # fcnName()

  fcnBody <- function(fcn) {
    paste(deparse(fcn), collapse="\n");
  } # fcnBody()

  # "Each [...] function evaluation has a tuple,
  #    (call, function definition, evaluation environment),
  #  coupled to it, which can be retrieved via sys.{call,function,frame}()."
  # Source: help("sys.parent", package="base")

  calls <- sys.calls();
  parents <- sys.parents();
  functions <- sys.functions(parents);

  stackTrace <- list();
  for (kk in seq_along(calls)) {
    call <- calls[[kk]];
    fcn <- functions[[kk]];
    name <- fcnName(call);
    body <- fcnBody(fcn);
    envir <- environment(fcn);
    envirName <- environmentName(envir);
    pathname <- fcnPathname(call);
    trace <- list(
      call=call,
      name=name,
      body=body,
      envir=envirName,
      pathname=pathname
    );
    stackTrace[[kk]] <- trace;
  } # for (kk ...)
  rm(list=c("calls", "parents", "functions"));

  # The new class is Exception, but for convenience it should also
  # derive from 'try-error', which is used by try() etc.
  extend(Object(), c("Exception", "simpleError", "error", "condition", "try-error"),
    .msg        = paste(..., sep=sep, collapse=collapse),
    .when       = Sys.time(),
    .stackTrace = stackTrace
  )
})




###########################################################################/**
# @RdocMethod as.character
#
# \title{Gets a character string representing of the Exception}
#
# \description{
#  @get "title".
#  By default the format is: "[\{POSIX date string\}] \{class name\}: \{msg\}".
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
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("as.character", "Exception", function(x, ...) {
  # To please R CMD check
  this <- x;

  paste("[", getWhen(this), "] ", class(this)[1], ": ", getMessage(this), sep="");
})



###########################################################################/**
# @RdocMethod print
#
# \title{Prints the Exception}
#
# \description{
#  @get "title". By default the \code{as.character()} representation plus
#  the stack trace is printed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{Returns nothing.}
#
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# @author
#
# \seealso{
#   @seemethod "as.character".
#   @seemethod "getStackTrace".
#   @seemethod "printStackTrace".
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("print", "Exception", function(x, ...) {
  cat(getStackTraceString(x, ...));
})




###########################################################################/**
# @RdocMethod getWhen
#
# \title{Gets the time when the Exception was created}
#
# \description{
#  Gets the time, as a POSIX object, when the Exception was created.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a POSIX time object.
# }
#
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("getWhen", "Exception", function(this, ...) {
  this$.when;
})





###########################################################################/**
# @RdocMethod getMessage
#
# @title "Gets the message of the Exception"
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
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("getMessage", "Exception", function(this, ...) {
  this$.msg;
})




###########################################################################/**
# @RdocMethod throw
#
# \title{Throws an Exception that can be caught}
#
# \description{
#  Throws an Exception that can be caught by \code{tryCatch()}.
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
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# @author
#
# \seealso{
#   @seeclass
#   See also \code{\link[base:conditions]{tryCatch}()}.
#   This method overrides (and is fully backward compatible with) the one
#   defined in the \pkg{R.methodsS3} package.
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("throw", "Exception", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Record this Exception
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Exception$.lastException <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Signal the exception as a condition
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  signalCondition(this);

  # Question: Are there ways to catch the above signals, and then via
  # some revoking mechanism continue below?!? /HB 2012-03-05


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If not caught by any handlers, output message containing the stack trace
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Output an error message containing the stacktrace
  msg <- getStackTraceString(this, ...);
  cat(msg, file=stderr());


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # ...and *abort*
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 1: Abort the current evaluation, but unfortunately abort()
  # is not guaranteed to not be "caught", cf. help("abort").
  # abort();

  # Alt 2: An alternative is to call stop() again, which will resignal
  # a condition and then abort.  The resignalled condition should not
  # really be caught by anything, because if so, it would have been
  # caught by the above signalling.  This is based on the assumption
  # that it is not possible to continue after the above signal,
  # iff it is caught. /HB 2012-03-05
  cond <- simpleCondition(getMessage(this));
  class(cond) <- "condition";
  stop(cond);
}, overwrite=TRUE, conflict="quiet")



###########################################################################/**
# @RdocMethod getLastException
#
# \title{Static method to get the last Exception thrown}
#
# \description{
#  Static method to get the last Exception instanciated.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "Exception" object.
# }
#
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# \seealso{
#   @seeclass
#   See also \code{\link[base:conditions]{tryCatch}()}.
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("getLastException", "Exception", function(this, ...) {
  Exception$.lastException;
}, static=TRUE);




###########################################################################/**
# @RdocMethod getStackTrace
#
# \title{Gets the stack trace saved when the exception was created}
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
#   Returns a @list containing the stack trace.
# }
#
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# \seealso{
#   @seemethod "printStackTrace".
#   \code{\link[utils:debugger]{dump.frames}()}.
#   \code{\link[base:conditions]{tryCatch}()}.
#   @seeclass
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("getStackTrace", "Exception", function(this, cleanup=getOption("R.oo::Exception/getStackTrace/args/cleanup", TRUE), ...) {
  stackTrace <- this$.stackTrace;
  names(stackTrace) <- seq_along(stackTrace);

  # Remove "uninformative" steps, e.g. tryCatch()
  if (cleanup) {
    # (a) Drop anything before doTryCatch()
    names <- sapply(stackTrace, FUN=function(x) x$name);
    idxs <- which(is.element(names, c("doTryCatch")));
    if (length(idxs) > 0) {
      idx <- idxs[length(idxs)];
      keep <- seq(from=idx+1L, to=length(stackTrace));
      stackTrace <- stackTrace[keep];
    }

    # (b) Drop anything after throw()
    names <- sapply(stackTrace, FUN=function(x) x$name);
    idxs <- which(is.element(names, "throw"));
    if (length(idxs) > 0) {
      idx <- idxs[1L];
      keep <- seq_len(idx-1L);
      stackTrace <- stackTrace[keep];
    }
  }

  stackTrace <- rev(stackTrace);

  stackTrace;
})



###########################################################################/**
# @RdocMethod getCalls
# @aliasmethod getCall
#
# \title{Gets the active calls saved when the exception was created}
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
#   Returns an unnamed @list with @language elements.
# }
#
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# \seealso{
#   @seemethod "getStackTrace".
#   \code{\link[utils:debugger]{sys.calls}()}.
#   @seeclass
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("getCalls", "Exception", function(this, ...) {
  stackTrace <- getStackTrace(this, ...);
  calls <- lapply(stackTrace, FUN=function(x) x$call);
  calls;
})

setMethodS3("getCall", "Exception", function(x, which=1L, ...) {
  # To please R CMD check (R >= 2.14.0)
  this <- x;
  calls <- getCalls(this, ...);
  if (length(calls) == 0) {
    return(NULL);
  }
  calls[[which]];
})




###########################################################################/**
# @RdocMethod getStackTraceString
#
# \title{Gets the stack trace as a string}
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
#   Returns a @character string.
# }
#
# \seealso{
#   @seemethod "getStackTrace".
#   @seeclass
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("getStackTraceString", "Exception", function(this, ..., details=TRUE) {
  head <- sprintf("%s\n", as.character(this));

  stackTrace <- getStackTrace(this, ...);

  if (length(stackTrace) == 0) {
    return(head);
  }

  calls <- sapply(stackTrace, FUN=function(trace) trace$call);
  res <- character(length=length(calls));
  for (kk in seq_along(calls)) {
    call <- calls[[kk]];
    rows <- deparse(call);
    if (details) {
      prefix1 <- sprintf("  at #%02d. ", as.integer(names(calls)[kk]));
    } else {
      prefix1 <- "  at ";
    }
    prefixT <- paste(rep(" ", times=nchar(prefix1)), collapse="");
    prefix <- rep(prefixT, times=length(rows));
    prefix[1] <- prefix1;
    rows <- sprintf("%s%s", prefix, rows);
    res[kk] <- paste(rows, collapse="\n");
  } # for (kk ...)

  if (details) {
    locs <- sapply(stackTrace, FUN=function(trace) {
      name <- trace$name;
      envir <- trace$envir;
      s <- sprintf("%s()", name);
      if (envir == "") {
        s <- sprintf("%s is local of the calling function", s);
      } else {
        s <- sprintf("%s is in environment '%s'", s, envir);
      }
      s;
    });

    res <- sprintf("%s\n          - %s", res, locs);

    pathnames <- sapply(stackTrace, FUN=function(trace) {
      trace$pathname;
    });
    pathnamesT <- sprintf("\n          - originating from '%s'", pathnames);
    pathnamesT[nchar(pathnames) == 0] <- "";
    res <- sprintf("%s%s", res, pathnamesT);

    res <- sprintf("%s\n", res);
  } # if (details)
  res <- paste(res, collapse="\n");

  if (details) {
    res <- sprintf("%s\n%s\n", head, res);
  } else {
    res <- sprintf("%s%s\n", head, res);
  }

  res;
}, private=TRUE)




###########################################################################/**
# @RdocMethod printStackTrace
#
# \title{Prints the stack trace saved when the exception was created}
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
#   Returns nothing.
# }
#
# \examples{\dontrun{For a complete example see help(Exception).}}
#
# \seealso{
#   @seemethod "getStackTrace".
#   \code{\link[base:conditions]{tryCatch}()}.
#   @seeclass
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("printStackTrace", "Exception", function(this, ...) {
  cat(getStackTraceString(this, ...));
})




############################################################################
# HISTORY:
# 2012-10-14
# o Now throw() for Exception outputs the error message both above and
#   below the stack trace, which is particularly useful when the stack
#   trace is long.
# 2012-09-14
# o ROBUSTNESS/BUG FIX: The Exception constructor could generate warning
#   'In if (regexpr("^function\\(", code) != -1) return("") : the
#   condition has length > 1 and only the first element will be used'
#   occuring in its local fcnName() function.  Now code no longer assumes
#   that 'code' is of length 1.
# 2012-09-10
# o Updated throw() for Exception to "abort" after signalling the condition
#   by calling stop() with an empty condition.  This is not perfect,
#   because it outputs an "Error:" message at the end, but it's better
#   than nothing.
# 2012-06-17
# o BUG FIX/GENERALIZATION: throw() for Exception would give an error on
#   R < 2.14.0, where no generic getCall() exists.  Now it works for
#   all versions of R.
# 2012-03-18
# o Now it is possible to set the default value of argument 'cleanup'
#   of getStackTrace() for Exception via an option.
# 2012-03-08
# o Now throw() for Exception utilizes abort().
# o Now Exception stores much more information about the stacktrace.
# 2012-03-07
# o Added getCalls() and getCall() for Exception.  Now Exception()
#   records the active calls, i.e. sys.calls().  getCalls() takes care
#   of the "cleanup" and "parsing".
# 2012-02-29
# o Now throw() of Exception utilizes new abort().
# o CLEANUP: Restructured the code of throw() for Exception.
# 2011-07-10
# o Changed first argument of getCall() to 'x', because that is what
#   the new getCall() method of 'stats' in R v2.14.0 uses.
# 2005-02-20
# o Updated broken link to tryCatch().
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2005-02-10
# o Moved showAndWait() from Exception to simpleError.
# 2004-10-18
# o Added more Rdoc comments.
# 2004-03-02
# o The Exception class now inherits from the simpleError, error and
#   condition classes.
# o The throw() method of Exception does now make use of the new (R v1.8.0)
#   signalCondition() method.
# 2003-12-16
# o Now throw() includes the complete stacktrace when generating an internal
#   error signal.
# 2003-04-13
# o Wrote Rdoc comments that were missing and updated some others.
# 2003-03-23
# o Added showAndAsk(), which will, if tcltk is installed, display a dialog
#   box with the error message. If tcltk is not installed, The message will
#   be printed on the command line and a prompt requesting the user to press
#   enter will be shown. showAndAsk() will give an error if run in a non-
#   interactive mode.
# 2003-01-19
# o Added stacktrace information to *each* Exception object. This is created
#   when the object is created.
# 2003-01-18
# o Replaced all occurences of getClass() with data.class(). Will change
#   the use of getClass() in the future to return a Class object.
# 2002-10-17
# o Made getLastException() a static method of Exception.
# o Created from previous ideas in R.oo.
############################################################################
