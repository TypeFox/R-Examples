
# general helper function to check arguments.
# if arg does not have proper class but can be coerced to it without changing
# the values, coerce the arg.
# otherwise, abort with appropriate error message:
# arg is not proper class
# arg is NULL and null.OK==FALSE
# length(arg) < min.length or length(arg) > max.length
check.arg <- function(arg, argname="argument", class=NULL,
                     null.OK=TRUE, min.length=1L, max.length=1L) {
  if (is.null(arg)) {
    if (null.OK) return(arg)
    stop(paste(argname, "cannot be NULL"))
  }

  # now arg is not null
  if ((!is.null(min.length) && length(arg) < min.length) ||
      (!is.null(max.length) && length(arg) > max.length)) {
    if (min.length == max.length)
      stop(paste(argname, "should be", class, "of length", min.length))
    stop(paste(argname, "should be ", class, "with length >=", min.length, "and length <=", max.length))
  }
  if (!is.null(class) &&
      !is(arg, class)) {
    y <- as(arg, class)
    if (sum(arg==y) == length(arg)) {
      return(y)
    }
  }
  
  if (!is.null(class) &&
      !is(arg, class))
    stop(paste(argname, "should be of type", class))
  return(arg)
}
