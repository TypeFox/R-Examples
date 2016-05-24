shiftDj <- function(ctmatrix){
shiftdj <- overallDiffCatj(ctmatrix) - quantityDj(ctmatrix) - exchangeDj(ctmatrix)
return(shiftdj)
}
