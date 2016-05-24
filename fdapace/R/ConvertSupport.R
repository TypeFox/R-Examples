# This function is used to convert the support of mu/phi/cov etc, to and from obsGrid and regGrid.
# If one wants to convert both mu and phi, it should be called one at a time.
# fromGrid and toGrid should be sorted.
# mu: any vector of a function
# phi: any matrix, each column containing a function to be interpolated.
# Cov: any matrix supported on fromGrid * fromGrid, to be interpolated to toGrid * toGrid. 

ConvertSupport <- function(fromGrid, toGrid, mu=NULL, Cov=NULL, phi=NULL) {

  # In case the range of toGrid is larger than fromGrid due to numeric error
  buff <- .Machine$double.eps * max(abs(fromGrid)) * 3
  if (abs(toGrid[1] - fromGrid[1]) < buff)
    toGrid[1] <- fromGrid[1]
  if (abs(toGrid[length(toGrid)] - fromGrid[length(fromGrid)]) < buff)
    toGrid[length(toGrid)] <- fromGrid[length(fromGrid)]
    
  if (!is.null(mu)) {# convert mu
    return(MapX1D(fromGrid, mu, toGrid))
  } else if (!is.null(Cov)) {
    gd <- pracma::meshgrid(toGrid) #pracma
    ret <- matrix(interp2lin(fromGrid, fromGrid, Cov, gd$X, gd$Y), nrow=length(toGrid))
    ret <- 0.5 * (ret + t(ret))                         # ensure that ret is symmetric
    return(ret)
  } else if (!is.null(phi)) {
    return(MapX1D(fromGrid, phi, toGrid))
  }

}

