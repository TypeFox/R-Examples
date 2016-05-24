fgl2vec <- function(vec){
  
  ## function to turn a state vector to the minimum numbering, e.g. s
  ## = c(2,3,1,1) to s = c(1,2,3,3).

  ## This ensures that when going from state -> label, the label is
  ## always the first label representing that state.
  
  y <- numeric( length(vec))
  xs <- unique(vec)
  for( i in 1:length(xs)){
    y <- replace(y, vec==xs[i], i)
  }
  return(y)
}
