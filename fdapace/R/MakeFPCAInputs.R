#' Format FPCA input
#'
#' Turn vector inputs to the list so they can be used in FPCA 
#' 
#' @param IDs  n-by-1 vector of subject IDs (Default: NULL)
#' @param tVec n-by-1 vector of measurement times
#' @param yVec n-by-1 vector of measurements from the variable of interest
#' @param na.rm logical indicating if NA should be omitted (Default: FALSE)
#' @return L list containing 3 lists each of length 'm', 'm' being the number of unique subject IDs
#' @export

MakeFPCAInputs <- function(IDs = NULL, tVec, yVec, na.rm=FALSE){

  if( !is.null(IDs) ){ 
    if (na.rm) {
      dat <- na.omit(data.frame(IDs, tVec, yVec))
      IDs <- dat[, 'IDs']
      tVec <- dat[, 'tVec']
      yVec <- dat[, 'yVec']
    }
    uniqueIDs <- unique(IDs) 
    Lt <- lapply( uniqueIDs, function(x) tVec[ which(IDs == x)])
    Ly <- lapply( uniqueIDs, function(x) yVec[ which(IDs == x)])
    Lid <- as.list(uniqueIDs)
  } else if ( is.matrix(yVec) && is.null(IDs) && is.vector(tVec) ){
    Ly <- lapply( seq_len(nrow(yVec)), function(i) yVec[i,])
    Lt <- rep( list(tVec), dim(yVec)[1] )
    Lid <- as.list( 1:dim(yVec)[1])
  }
  L <- list( Lid = Lid, Ly = Ly, Lt = Lt)
  return(L)

}

