# Function: families
# Argument: Multiple character strings.
# Description: This function is used mostly for the user's convenience. It simply creates a list of character strings and
# can be used in the specification of parameters for gatekeeping procedures.
#' @export
families = function(...) {

  args = list(...)

  nargs = length(args)

  if (nargs <= 0) stop("Families function: At least one family must be specified.")

  return(args)
  invisible(args)

}