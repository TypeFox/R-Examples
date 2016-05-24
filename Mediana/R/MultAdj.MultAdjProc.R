######################################################################################################################

# Function: MultAdj.MultAdjProc.
# Argument: MultAdjProc object
# Description: This function is used to create an object of class MultAdjProc.
#' @export
MultAdj.MultAdjProc = function(...) {


  multadj = lapply(list(...),function(x) {if(class(x)=="MultAdjProc") list(unclass(x)) else unclass(x)} )

  class(multadj) = "MultAdj"
  return(multadj)
  invisible(multadj)

}