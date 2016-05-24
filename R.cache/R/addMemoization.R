#########################################################################/**
# @RdocDefault addMemoization
#
# @title "Creates a copy of an existing function such that its results are memoized"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fcn}{A @function (or the name of a function) that should be
#     copied and have memoization added.}
#   \item{envir}{The @environment from where to look for the function.}
#   \item{...}{Additional arguments for controlling the memoization,
#     i.e. all arguments of @see "memoizedCall" that are not passed
#     to @see "base::do.call".}
# }
#
# \value{
#   Returns a @function.
# }
#
# \details{
#  The new function is setup such that the the memoized call is done
#  in the environment of the caller (the parent frame of the function).
#
#   If the @function returns @NULL, that particular function call is
#   \emph{not} memoized.
# }
#
# @author
#
# \seealso{
#  The returned function utilized @see "memoizedCall" internally.
# }
#
# @keyword "programming"
# @keyword "IO"
#*/#########################################################################
setMethodS3("addMemoization", "default", function(fcn, envir=parent.frame(), ...) {
  # Argument 'fcn':
  if (is.character(fcn)) {
    if (!exists(fcn, mode="function", envir=envir, inherits=TRUE)) {
      throw("Argument 'fcn' is not an existing function: ", fcn);
    }
    fcn <- get(fcn, mode="function", envir=envir, inherits=TRUE);
  }

  if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", mode(fcn));
  }

  # Already memoized?
  if (inherits(fcn, "MemoizedFunction")) {
    return(fcn)
  }

  # Record the argument specific to memoizedCall().
  memArgs <- list(...);

  res <- function(..., envir=parent.frame()) {
    args <- list(fcn, ..., envir=envir);
    args <- c(args, memArgs);
    do.call("memoizedCall", args=args);
  }
  class(res) <- c("MemoizedFunction", class(res))

  res
}) # addMemoization()


#######################################################################
# HISTORY:
# 2014-09-10
# o ROBUSTNESS: addMemoization() will no longer memoize an already
#   memoized function.
# 2011-02-14
# o Added addMemoization().
#######################################################################
