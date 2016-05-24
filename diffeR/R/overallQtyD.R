overallQtyD <- function(ctmatrix){
qtydj <- abs(apply(ctmatrix, 1, sum) - apply(ctmatrix, 2, sum))
overallqtyd <- sum(qtydj)/2
return(overallqtyd)
}
