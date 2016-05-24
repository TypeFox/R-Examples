quantityDj <- function(ctmatrix){
qtydj <- abs(apply(ctmatrix, 1, sum) - apply(ctmatrix, 2, sum))
return(qtydj)
}
