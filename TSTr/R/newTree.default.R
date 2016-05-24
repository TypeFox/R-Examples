#' @export
#' 
newTree.default <-
function(input) {
  arbol <- tstTree()
  contador <- 0
  if(!file.exists(input)[1]){
    for (word in input) {
      arbol <- addWord(arbol,as.character(word))
      contador <- contador + 1
      if (contador %% 10000 == 0) {
        print(contador)
      }
    }
  }
  else if (is.null(input)) {return()}
  else{
    content <- readLines(input, warn = FALSE)
    for (word in content) {
      arbol <- addWord(arbol,as.character(word))
      contador <- contador + 1
      if (contador %% 10000 == 0) {
        print(contador)
      }
    }
  }
  cat(paste("Tree created with",contador,"words and",length(arbol$ch),"nodes\n"))
  return(arbol)
}
