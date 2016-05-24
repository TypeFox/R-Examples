adj.measlev <-
function(x,threshold=.5){
  
  x@.Data <- lapply(x@.Data, adj.input, threshold=threshold)
  x
}
