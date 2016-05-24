# This function converts dense regular functional input list
# to a matrix for easier dense case implementation
##########################################################################
# Input:  - y: list of n dense regular observed p-dim functional objects
##########################################################################
# Output: - ymat: n by p matrix containing all functional data
##########################################################################

List2Mat <- function(y,t){
  n = length(y)
  obsGrid = sort(unique(unlist(t)))
  ymat = matrix( rep(NA, n * length(obsGrid)), nrow = n, byrow = TRUE)
  
  for (i in 1:n){
    ymat[i, is.element(obsGrid, t[[i]])] = y[[i]]   
  }
  return(ymat)
}
