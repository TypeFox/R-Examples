exchangeDij <- function(ctmatrix){
exchDij <- 2*pmin(ctmatrix, t(ctmatrix))
exchDij.lt <- lower.tri(exchDij)*exchDij
return(exchDij.lt)
}
