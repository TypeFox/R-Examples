#' Sparsify densely observed functional data
#'
#' Given a matrix of densely observed functional data, make a sparsified sample.
#' 
#' @param samp A matrix of densely observed functional data, with each row containing one sample.
#' @param pts A vector of grid points corresponding to the columns of \code{samp}.
#' @param sparsity A vector of integers. The number of observation per sample is chosen to be one of the elements in sparsity with equal chance.
#' @param aggressive Sparsify in an "aggressive" manner making sure that near-by readings are excluded.
#' @param fragment Sparsify the observations into fragments, which are (almost) uniformly distributed in the time domain. Default to FALSE as not fragmenting. Otherwise a positive number specifying the approximate length of each fragment.
#'
#' @return A list of length 2, containing the following fields:
#' \item{tList}{A list of observation time points for each sample.}
#' \item{yList}{A list of values for each sample, corresponding to the time points.}
#' @export
Sparsify <- function(samp, pts, sparsity, aggressive = FALSE, fragment=FALSE) {

    if (length(sparsity) == 1)
      sparsity <- c(sparsity, sparsity) # avoid scaler case

    if (aggressive && fragment)
      stop('Specify one of `aggressive` or `fragment` only')
      
    if (aggressive) {
      indEach <- lapply(1:nrow(samp), function(x)
        remotesampling(ncol(samp), sparsity) )
    } else if (fragment != FALSE) {
      nptsEach <- fragment / mean(diff(pts))
      indEach <- lapply(1:nrow(samp), function(x) {
        ranPts <- range(pts)
        mid <- runif(1, ranPts[1], ranPts[2])
        usePts <- which(pts >= mid - 1/2 * diff(ranPts) * fragment & 
                        pts <= mid + 1/2 * diff(ranPts) * fragment)
        nSampPts <- sample(sparsity, 1)
        if (nSampPts >= length(usePts)) usePts else sort(sample(usePts, nSampPts))
      })
    } else {
      indEach <- lapply(1:nrow(samp), function(x) 
        sort(sample(ncol(samp), sample(sparsity, 1))))
    }
    
    tList <- lapply(indEach, function(x) pts[x])
    yList <- lapply(1:length(indEach), function(x) {
        ind <- indEach[[x]]
        y <- samp[x, ind]
        return(y)
    })
   
    return(list(tList=tList, yList=yList))
}



remotesampling <- function(N,s){
  onesamp = sort( sample( N, sample( s, 1))) 
  threshold = (1/length(onesamp))^(1.5) * N 
  while( min(diff(onesamp)) < threshold ){
    onesamp = sort( sample( N, length(onesamp)))
  }
  return(onesamp)
}



