#' Create lagged RasterStacks
#' 
#' @description
#' The function is used to produce two lagged RasterStacks. The second is cut
#' from the beginning, the first from the tail to ensure equal output lengths
#' (provided that input lengths were equal).
#' 
#' @param x a RasterStack (to be cut from tail)
#' @param y a RasterStack (to be cut from beginning)
#' @param lag the desired lag (in the native frequency of the RasterStack)
#' @param freq the frequency of the RasterStacks
#' @param ... currently not used
#' 
#' @return
#' a list with the two RasterStacks lagged by \code{lag}
#' 
#' @examples
#' data(pacificSST)
#' data(australiaGPCP)
#' 
#' # lag GPCP by 4 months
#' lagged <- lagalize(pacificSST, australiaGPCP, lag = 4, freq = 12)
#' lagged[[1]][[1]] #check names to see date of layer
#' lagged[[2]][[1]] #check names to see date of layer
#' 
#' @export lagalize
lagalize <- function(x, 
                     y, 
                     lag = NULL,
                     freq = 12,
                     ...) {
  
  # Return list of unmodified RasterStacks if lag == NULL
  if (is.null(lag)) {
    return(list(x, y))
  } else {
    rest <- freq - lag
    
    # Lagalize predictor stack
    x.lag <- cutStack(x = x, tail = TRUE, n = lag)
    x.lag.adj <- x.lag[[1:(raster::nlayers(x.lag) - rest)]]
    
    # Lagalize response stack
    y.lag <- cutStack(x = y, tail = FALSE, n = lag)
    y.lag.adj <- y.lag[[1:(raster::nlayers(y.lag) - rest)]]
    
    # Return list of lagalized stacks
    return(list(x.lag, y.lag))
  }
  
}