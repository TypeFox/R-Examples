#' @export
#' 
SDkeeper.default <- function(input, maxdist, useTST = FALSE){
  if(useTST==F){
    crearFullDT(input, maxdist)
  }
  else{
    crearFullTree(input, maxdist)
  }
}