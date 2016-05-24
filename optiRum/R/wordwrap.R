#' Produce a string with one word per line
#' Designed for splitting strings to fit better on axis on charts
#' 
#' @param x string
#' @param ... Allows additional parameters to be passed to gsub
#' 
#' @keywords ggplot
#' @family helper  
#' 
#' @export
#' 
#' @examples
#' library('ggplot2')
#' ggplot(data.frame(x=1:10,y=1:10),aes(x,y))+theme_optimum()+geom_line()
#' 
wordwrap <- function(x, ...) {
    gsub(pattern = "\\s", replacement = "\n", x, ...)
} 
