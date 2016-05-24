spa <- function(Y,X,C){
  both <- suppressWarnings(ma(cbind(Y=Y,X=X,C=C),partition=list(c(1),c(2,3)))$A)
  one <- suppressWarnings(ma(cbind(Y=Y,C=C),partition=list(c(1),c(2)))$A)
  return (max(both-one,0))
}
