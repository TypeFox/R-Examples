#' Un-rowname
#' 
#' Strip rownames from an object (stolen from plyr).
#' 
#' @author Hadley Wickham
#' @keywords keywords
#' 
#' @param x A data.frame.
#'   
#' @return A data.frame with rownames removed.
#' 
#' @examples
#' unrowname(mtcars)
#'   
#' @export

unrowname <- function(x) {
    rownames(x) <- NULL
    x
}
