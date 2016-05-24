nearestToArchetypes <- function(indivs,numArch,mdras){
 as.numeric(which.min(mdras[indivs,]) - (numArch))
}
