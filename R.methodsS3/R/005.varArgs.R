hasVarArgs <- function(...) UseMethod("hasVarArgs");
export(hasVarArgs) <- TRUE;

hasVarArgs.function <- function(fcn, ...) {
  if (!is.function(fcn))
    stop("Argument 'fcn' must be a function: ", mode(fcn));

  # Get the current formals
  args <- formals(fcn);

  is.element("...", names(args));
} # hasVarArgs()
S3class(hasVarArgs.function) <- "function";
export(hasVarArgs.function) <- FALSE;


appendVarArgs <- function(...) UseMethod("appendVarArgs");
export(appendVarArgs) <- TRUE;

appendVarArgs.function <- function(fcn, ...) {
  if (hasVarArgs(fcn))
    return(fcn);

  # Get the current formals
  args <- formals(fcn);
  # Add '...'
  args <- c(args, formals(function(...) {}));
  # Set new formals
  formals(fcn) <- args;

  fcn;
} # appendVarArgs()
S3class(appendVarArgs.function) <- "function";
export(appendVarArgs.function) <- FALSE;


############################################################################
# HISTORY:
# 2005-02-15
# o Created.
############################################################################
