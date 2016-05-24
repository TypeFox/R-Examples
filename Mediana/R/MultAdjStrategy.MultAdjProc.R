######################################################################################################################

# Function: MultAdjStrategy.MultAdjProc.
# Argument: MultAdjProc object
# Description: This function is used to create an object of class MultAdjStrategy.
#' @export
MultAdjStrategy.MultAdjProc = function(...) {

  multadjstrat = lapply(list(...),unclass)

  class(multadjstrat) = "MultAdjStrategy"
  return(multadjstrat)
  invisible(multadjstrat)

}