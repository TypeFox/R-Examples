#' @export
#' 
dimTree.tstTree <-
function(tree){
  return(c(sum(tree$flag),length(tree$ch)))
}
