print.diagram <-
function(x, ...) {
  n <- nrow(x)
  print(x[seq_len(n), ])  
}