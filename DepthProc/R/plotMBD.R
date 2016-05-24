#' Functional boxplot based on Modified Band Depth
#' 
#' @param x object of class DepthMBD
#' @param band_lim limits for bands
#' @param colors color set
#' @param add_lines Logical. If TRUE, all lines from data will be plotted.
#' @param lwd lines width
#' @param alpha transparency for band area.
#' @param \dots other arguments passed to matplot
#'
#' @examples
#' x = depthMBD(matrix(rnorm(200), ncol = 10))
#' plot(x)
#'
setMethod("plot", signature = c(x = "DepthMBD"), function(x, band_lim = c(0,0.5), colors = NULL, add_lines = TRUE, lwd = 1, alpha = 0.5, ...)
{
  tu = t(x@u)
  
  matplot(tu, type = "n", col = "gray90", lty = 1, ...)
  
  band_lim = quantile(x, band_lim)
  if(is.null(colors))
  {
    colors = rev(gray.colors(length(band_lim), start = 0.4))
  }
  
  col_lines = colors
  if(alpha != 1) 
  {
    colors = .addAlpha(colors, alpha)
  }
  
  for(i in seq_along(band_lim))
  {
    
    tmp_tu = tu[,x >= band_lim[i]]
    range = t(apply(tmp_tu, 1, range))
    range[,2] = rev(range[,2])
    xx = c(1:nrow(range), nrow(range):1)
    polygon(xx, range, border = NA, col = colors[i])
    
  }
  
  if(add_lines)
  {
    col_lines
    col = rep(0,ncol(tu))
    
    for(i in seq_along(band_lim))
    {
      col = col + as.numeric(x >= band_lim[i])
    }
    
    sapply(1:ncol(tu), function(i) lines(tu[,i], col = col_lines[col[i]]))
    lines(tu[,which.max(x)], lwd = lwd*2)
  }
  
  
  
  

  
})


.addAlpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, function(x) 
    rgb(x[1], x[2], x[3], alpha=alpha))  
}



