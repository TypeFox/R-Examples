#' @export
#' 
searchWord.tstTree <-
function(tree, string){
  result <- 0
  node <- 1
  i <- 1
  while (!is.na(node)) {
    if (str_sub(string,i,i) < tree$ch[node]) {
      node <- tree$L[node]
    }
    else if (str_sub(string,i,i) > tree$ch[node]) {
      node <- tree$R[node]
    }
    else{
      i <- i + 1
      if (i - 1 == nchar(string)) {
        result <- tree$flag[node]
      }
      else{
        node <- tree$C[node]
      }
    }
  }
  as.logical(result)
}
