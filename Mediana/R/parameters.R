# Function: parameters
# Argument: Multiple character strings.
# Description: This function is used mostly for the user's convenience. It simply creates a list of character strings and
# can be used in cases where multiple parameters need to be specified.
#' @export
parameters = function(...) {

  args = list(...)

  nargs = length(args)

  if (nargs <= 0) stop("Parameters function: At least one parameter must be specified.")

  return(args)
  invisible(args)

}