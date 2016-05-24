overallDiff <- function(ctmatrix){
overalldiff <- sum(ctmatrix) - sum(diag(ctmatrix))
return(overalldiff)
}