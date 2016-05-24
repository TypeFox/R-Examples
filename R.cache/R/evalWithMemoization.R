#########################################################################/**
# @RdocFunction evalWithMemoization
#
# @title "Evaluates an R expression with memoization"
#
# \description{
#  @get "title" such that the same objects are assigned to the
#  current environment and the same result is returned, if any.
# }
#
# @synopsis
#
# \arguments{
#   \item{expr}{The @expression to be evaluated.}
#   \item{key}{Additional objects to uniquely identify the evaluation.}
#   \item{...}{Additional arguments passed to @see "loadCache" 
#     and @see "saveCache".}
#   \item{envir}{The @environment in which the expression should
#     be evaluated.}
#   \item{force}{If @TRUE, existing cached results are ignored.}
# }
#
# \value{
#   Returns the value of the evaluated \code{expr} @expression, if any.
# }
#
# @examples "../incl/evalWithMemoization.Rex"
#
# @author
#
# \seealso{
#  Internally, @see "base::eval" is used to evaluate the expression.
# }
#
# @keyword "programming"
# @keyword "IO"
#*/#########################################################################  
evalWithMemoization <- function(expr, key=NULL, ..., envir=parent.frame(), force=FALSE) {
  expr <- substitute(expr);

  # Setup a unique list of keys
  key <- c(list(expr=expr), key);

  # Look for cached results
  resList <- loadCache(key=key, ...);
  if (!force && !is.null(resList)) {
    # Attach all objects memoized during the evaluation
    attachLocally(resList$envir, envir=envir);
    # Return the results of the memoized evaluation
    return(resList$result);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Evaluate expression
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Evaluate the expression in a temporary environment, so that
  # we memoize all objects created along with the results.
  env <- new.env(parent=envir);
  res <- eval(expr, envir=env);

  # NOTE: For some unknown reason does attachLocally() set 
  # the fields inside 'env' to NULL.  /HB 2011-04-02
  fields <- ls(envir=env, all.names=TRUE);
  for (field in fields) {
    assign(field, get(field, envir=env), envir=envir);
  }

  # Cache results
  resList <- list(envir=env, results=res);
  saveCache(resList, key=key, ...);

  res;
} # evalWithMemoization()


#######################################################################
# HISTORY:
# 2011-04-01
# o Added evalWithMemoization().
#######################################################################
