#' @title ggplot2 tools
#'
#' @description A few handy tools for working with ggplot2.
#'
#' @details
#' The \code{ggplot2_extract_legend} function returns a list with the first
#' element being the legend and the second the original plot with the legend
#' omitted.
#'
#' @param x a ggplot
#'
#' @return a list with each element a ggplot
#' 
#' @examples
#' # a simple plot
#' my_plot <- 
#'   ggplot2::ggplot(mtcars) + 
#'   ggplot2::aes(x = wt, y = mpg, color = wt, shape = factor(cyl)) + 
#'   ggplot2::geom_point()
#' 
#' my_plot
#' 
#' # extract the legend.  the return object is a list with two elements, the first
#' # element is the legend, the second is the original plot sans legend.
#' temp <- ggplot2_extract_legend(my_plot)
#' 
#' # view just the legend.  This can be done via a call to the object or using
#' # plot or print.
#' temp
#' plot(temp[[1]])
#' 
#' # the original plot without the legened
#' plot(temp[[2]])
#'
#' @export   
#' @rdname ggplot2_tools

ggplot2_extract_legend <- function(x) { 
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  lgnd <- tmp$grobs[[leg]] 

  out <- list(legend = lgnd,
              plot   = x + ggplot2::theme(legend.position = "none"))
  class(out) <- c("qwraps2_ggplot2_extracted_legend", class(out))
  out
}

#' @export
print.qwraps2_ggplot2_extracted_legend <- function(x, ...) {
  graphics::plot(x[[1]])
}





