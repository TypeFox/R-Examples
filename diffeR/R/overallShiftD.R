overallShiftD <- function(ctmatrix){
Shift <- shiftDj(ctmatrix)
overallshfd <- sum(Shift)/2
return(overallshfd)
}
