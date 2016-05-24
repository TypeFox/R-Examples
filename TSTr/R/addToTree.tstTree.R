#' @export
#' 
addToTree.tstTree <-
function(tree, input) {
  contador <- 0
  if (!file.exists(input)[1]){
    for (word in input) {
      tree <- addWord(tree, as.character(word))
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
      tree <- addWord(tree,as.character(word))
      contador <- contador + 1
      if (contador %% 10000 == 0) {
        print(contador)
      }
    }
  }
  cat(paste(" ",contador,"words added and the total nodes are",length(tree$ch),"\n"))
  return(tree)
}
