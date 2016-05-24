######################################################################################################################

# Function: MultAdjProc.
# Argument: ....
# Description: This function is used to call the corresponding function according to the class of the argument.
#' @export
MultAdjProc = function(proc, par = NULL, tests = NULL) {

  # Error checks
  if (!is.na(proc) & !is.character(proc)) stop("MultAdj: multiplicity adjustment procedure must be character.")
  if (!is.null(par) & !is.list(par)) stop("MultAdj: par must be wrapped in a list.")
  if (!is.null(tests) & !is.list(tests)) stop("MultAdj: tests must be wrapped in a list.")
  if (any(lapply(tests, is.character) == FALSE)) stop("MultAdj: tests must be character.")


  mult.adjust = list(proc = proc, par = par, tests = tests)

  class(mult.adjust) = "MultAdjProc"
  return(mult.adjust)
  invisible(mult.adjust)
}