# Function: statistics
# Argument: Multiple character strings.
# Description: This function is used mostly for the user's convenience. It simply creates a list of character strings and
# can be used in cases where multiple statistics need to be specified.
#' @export
statistics = function(...) {

  args = list(...)

  nargs = length(args)

  if (nargs <= 0) stop("Statistics function: At least one test must be specified.")

  return(args)
  invisible(args)

}