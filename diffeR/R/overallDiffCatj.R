overallDiffCatj <- function(ctmatrix){
overalldiffj <- apply(ctmatrix, 1, sum) + apply(ctmatrix, 2, sum) - 2*diag(ctmatrix)
return(overalldiffj)
}