#' mlxR wrapper for ggplot
#' 
#' mlxR wrapper around \code{\link[ggplot2]{ggplot}} with a custom theme
#' 
#' @param \dots parameters passed to \code{\link[ggplot2]{ggplot}}
#' 
#' @return see \code{\link[ggplot2]{ggplot}}
#' @importFrom ggplot2 ggplot theme_bw theme element_rect
#' @importFrom grDevices rgb
#' @export         
ggplotmlx <- function(...) {
  ggplot(...) + theme_bw() + theme(plot.background = element_rect(fill=rgb(1,1,1))) 
}

