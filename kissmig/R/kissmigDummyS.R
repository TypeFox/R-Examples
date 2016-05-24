#' Get suitability maps of Europe 
#'
#' \command{kissmigDummyS} is a support function to generate suitability maps of Europe for example code.
#' @usage
#' kissmigDummyS(mean, sd, download=FALSE)
#' @param mean temperature mean (degree celsius) of the suitability distribution
#' @param sd temperature standard deviation (degree celsius) of the  suitability distribution
#' @param download if TRUE, required climate data are downloaded from www.worldclim.org
#'
#' @details
#' \command{kissmigDummyS} returns a suitability map of Europe based on mean annual temperature. It uses
#' data of worldclim and calculates suitability as a normal distribution defined by \command{mean} and \command{sd}
#' of mean annual temperature. The density function is linarely rescaled to a maximum of 1, the maximum suitability
#' used in \command{kissmig}. Set \command{download=TRUE} to download the required climate data when running the
#' function the first time.
#' @seealso \code{\link{kissmig}}
#' @export kissmigDummyS  
#' @references \url{www.worldclim.org}  

kissmigDummyS <- function(mean, sd, download=FALSE) {
  wcl <- getData('worldclim', var='bio', res=5, download=download) # download worldclim data
  ans <- crop(raster(wcl, layer=1), extent(-12,38,28,72))          # mean annual temperature for Europe
  if (((mean*10)<minValue(ans)) | ((mean*10)>maxValue(ans))) {
    warning("mean outside temperature range")
  }
  mean <- mean*10; sd <- sd*10
  values(ans) <- dnorm(values(ans), mean, sd)/dnorm(mean, mean, sd)
  return(ans)
}
