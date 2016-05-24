# Function: tests
# Argument: Multiple character strings.
# Description: This function is used mostly for the user's convenience. It simply creates a list of character strings and
# can be used in cases where multiple tests need to be specified.
#' @export
tests = function(...) {

  args = list(...)

  nargs = length(args)

  if (nargs <= 0) stop("Tests function: At least one test must be specified.")

  return(args)
  invisible(args)

}