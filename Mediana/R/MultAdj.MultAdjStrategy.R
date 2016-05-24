######################################################################################################################

# Function: MultAdj.MultAdjStrategy.
# Argument: MultAdjStrategy object
# Description: This function is used to create an object of class MultAdjStrategy.
#' @export
MultAdj.MultAdjStrategy = function(...) {

  multadj = lapply(list(...),function(x) {if(class(x)=="MultAdjProc") list(unclass(x)) else unclass(x)} )

  class(multadj) = "MultAdj"
  return(multadj)
  invisible(multadj)

}