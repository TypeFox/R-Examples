#' Build a character vector or list with adaptive string formatting
#'
#' The \code{rprintf} function checks the given character vector or
#' list and applies appropriate formatters that transform it from
#' generic patterns to specific texts with variables and indices
#' as placeholders replaced by a given set of values in correct
#' formats.
#'
#' @param x The character vector or list to be transformed
#' @param ... The arguments that specify the set of values to be
#'   placed
#' @export
#' @examples
#' \dontrun{
#'#' # Format a single-entry character vector with sprintf mechanism
#' rprintf('Hello, %s','world')
#' rprintf('%s (%d years old)','Ken',24)
#' rprintf('He is %d but has a height of %.1fcm',18,190)
#'
#' # Format a single-entry character vector with variable mechanism
#' rprintf('Hello, $name', name='world')
#' rprintf('$name ($age years old)',name='Ken',age=24)
#' rprintf('He is $age but has a height of $height:.2fcm',age=18,height=190)
#' rprintf('$a, $b:.1f, $c:+.2f, $b, $a:.0f',a=1.56,b=2.34,c=3.78)
#'
#' # Format a single-entry character vector with numbering mechanism
#' rprintf('Hello, {1}', 'world')
#' rprintf('{1} ({2} years old)','Ken',24)
#' rprintf('He is {1} but has a height of {2:.2f}cm',18,190)
#' rprintf('{1}, {2:.1f}, {3:+.2f}, {2}, {1:.0f}',1.56,2.34,3.78)
#' rprintf('{2},{1}','x','y')
#'
#' # This function also works for character vectors and lists.
#' rprintf(c('%s:%d','$name:$age','{1}:{2}'),name='Ken',age=24)
#' rprintf(c(a='%s:%d',b='$name:$age',c='{1}:{2}'),name='Ken',age=24)
#' rprintf(list('%s:%d','$name:$age','{1}:{2}'),name='Ken',age=24)
#' rprintf(list(a='%s:%d',b='$name:$age',c='{1}:{2}'),name='Ken',age=24)
#'
#' # It also works with list argument for named variables.
#' p <- list(name='Ken',age=24)
#' rprintf('name: $name, age: $age',p)
#' rprintf('name: {1}, age: {2}',p)
#'
#' Note that when the list of arguments are given names,
#' the variable names in format string should be modified.
#' rprintf('name: $arg.name, age: $arg.age', arg = p)
#' }
rprintf <- function(x, ...) {
  matches <- do.call(cbind, lapply(patterns, function(p) grepl(p, x, perl = TRUE)))
  funs.id <- apply(matches, 1L, function(row) which(row)[1L])
  funs <- names(patterns)[funs.id]
  result <- Map(rprintf.match, x, funs, list(list(...)))
  if (is.list(x)) 
    result else setnames(unlist(result, use.names = FALSE), names(x))
} 
