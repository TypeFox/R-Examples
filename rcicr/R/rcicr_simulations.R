#' Simulate pixel intensity range for noise 
#' 
#' @export
#' @import dplyr
#' @import matlab
#' @importFrom utils data
#' @importFrom stats runif
#' @importFrom graphics boxplot
#' @param nrep Number of replications
#' @param img_size Size of noise pattern in pixels (one value equal for width and height)
#' @return Matrix with range of noise intensities for each replication 
simulateNoiseIntensities <- function(nrep=1000, img_size=512) {
  
  results <- matlab::zeros(nrep, 2)
  s <- generateNoisePattern(img_size=512)
  
  pb <- dplyr::progress_estimated(length(unique(data[,by])))
  for (i in 1:nrep) {
    pb$tick()$print()
    
    params <- (runif(4096) * 2) - 1
    
    noise <- generateNoiseImage(params, s) 
    results[i,] <- range(noise)
  }
  pb$stop()
  boxplot(results)
  return(results)
}
