#' Shorten a RasterStack
#' 
#' @description
#' The function cuts a specified number of layers off a RrasterStack in 
#' order to create lagged RasterStacks.
#' 
#' @param x a RasterStack
#' @param tail logical. If \code{TRUE} the layers will be taken off
#' the end of the stack. If \code{FALSE} layers will be taken off
#' the beginning.
#' @param n the number of layers to take away.
#' 
#' @return a RasterStack shortened by \code{n} layers either from the 
#' beginning or the end, depending on the specification of \code{tail}
#' 
#' @examples
#' data(australiaGPCP)
#' 
#' # 6 layers from the beginning
#' cutStack(australiaGPCP, tail = FALSE, n = 6)
#' # 8 layers from the end
#' cutStack(australiaGPCP, tail = TRUE, n = 8)
#' 
#' @export cutStack
cutStack <- function(x, 
                     tail = TRUE,
                     n = NULL) {
  
  # Return unmodified RasterStack if n == NULL
  if (is.null(n)) {
    return(x)
  } else {
    # Supplied RasterStack is predictor:
    if (tail) {
      return(x[[1:(raster::nlayers(x)-n)]])
    # Supplied RasterStack is response:  
    } else {
      return(x[[(n+1):raster::nlayers(x)]])
    }
  }
  
}