# Function: samples
# Argument: Multiple character strings.
# Description: This function is used mostly for user's convenience. It simply creates a list of character strings and
# can be used in cases where multiple samples need to be specified.
#' @export
samples = function(...) {

  args = list(...)

  nargs = length(args)

  if (nargs <= 0) stop("Samples function: At least one sample must be specified.")

  return(args)
  invisible(args)

}