##create a trend
trend <- cbind(1:5,sin(1:5))
##an index of locations
idx <- c(rep(1:3,3),1:2,2:3)
##a list of time points for each location/observation
T <- c(rep(1:3,each=3),4,4,5,5)

##expand the F matrix to match the locations/times in idx/T.
F <- trend[T,]

##compute the expanded matrix
expandF(F, idx)

##compute the expanded matrix, assuming additional locations
expandF(F, idx, 5)

##or as a full matrix
expandF(F, idx, 5, sparse=FALSE)

\dontshow{
  if( abs(max(expandF(F, idx)-expandF(F, idx, sparse=FALSE)))>1e-13 ){
    stop("Error in 'expandF', full not equal")
  }
  if( abs(max(expandF(F, idx, 5)-expandF(F, idx, 5, sparse=FALSE)))>1e-13 ){
    stop("Error in 'expandF', expanded full not equal")
  }
}
