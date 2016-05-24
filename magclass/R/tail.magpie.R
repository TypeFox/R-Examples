tail.magpie <- function(x, n1=3L, n2=6L, n3=2L, ...) {
  if(dim(x)[1]<n1) n1 <- dim(x)[1]
  if(dim(x)[2]<n2) n2 <- dim(x)[2]
  if(dim(x)[3]<n3) n3 <- dim(x)[3]
  return(x[dim(x)[1]+1-(n1:1),dim(x)[2]+1-(n2:1),dim(x)[3]+1-(n3:1)])  
}