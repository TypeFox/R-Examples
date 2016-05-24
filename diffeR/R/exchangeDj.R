exchangeDj <- function(ctmatrix){
exchDij <- 2*pmin(ctmatrix, t(ctmatrix))
exchDij.lt <- lower.tri(exchDij)*exchDij
exchDj <- apply(exchDij.lt, 1, sum) + apply(exchDij.lt, 2, sum)
return(exchDj)
}
