charEdgelist <-
function(edgelist, vertexNames) {
  char <- matrix('0', ncol = ncol(edgelist), nrow = nrow(edgelist))

  for (col in 1:ncol(edgelist)){
    for (row in 1:nrow(edgelist)) {
      char[row, col] <- vertexNames[edgelist[row, col]]
    }
  }
  return(char)
}
