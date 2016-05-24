#' Create seasonal anomalies
#' 
#' @description
#' The function calculates anomalies of a RasterStack by supplying a 
#' suitable seasonal window. E. g. to create monthly anomalies of a 
#' raster stack of 12 layers per year, use \code{cycle.window = 12}.
#' 
#' @param x a RasterStack
#' @param cycle.window the window for the creation of the anomalies
#' @param use.cpp logical. Determines whether to use \strong{Rcpp} 
#' functionality, defaults to \code{TRUE}.
#' @param ... currently not used
#' 
#' @return a deseasoned RasterStack
#' 
#' @seealso
#' \code{\link{anomalize}}, \code{\link{denoise}}
#' 
#' @export deseason
#' 
#' @examples 
#' data("australiaGPCP")
#' 
#' aus_dsn <- deseason(australiaGPCP, 12)
#' 
#' opar <- par(mfrow = c(1,2))
#' plot(australiaGPCP[[1]], main = "original")
#' plot(aus_dsn[[1]], main = "deseasoned")
#' par(opar)
deseason <- function(x, 
                     cycle.window = 12,
                     use.cpp = FALSE,
                     ...) {
  
  if (use.cpp) {
    ## raster to matrix
    mat <- as.matrix(x)
    
    ## deseasoning
    mat_mv <- monthlyMeansC(mat, cycle.window)
    x_mv <- x[[1:cycle.window]]
    x_mv <- setValues(x_mv, values = mat_mv)
  } else {
    # Calculate layer averages based on supplied seasonal window
    x_mv <- raster::stack(rep(lapply(1:cycle.window, function(i) {
      raster::calc(x[[seq(i, raster::nlayers(x), cycle.window)]], fun = mean)
    }), raster::nlayers(x) / cycle.window))
  }
  
  # Subtract monthly averages from actually measured values
  x_dsn <- x - x_mv
  
  # Return output
  return(x_dsn)
}
