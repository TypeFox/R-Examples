#' Retrieve detail results of univariate regression models
#'


#' @param x A reg object
#' @export
#' @seealso \code{\link{reg}}
#' @examples
#' data(diabetes)
#' head(diabetes)
#'
#' reg_glm<-reg(data = diabetes, y = 5, factor = c(1, 3, 4), model = 'glm')
#' detail(reg_glm)


detail.reg <- function(x) {
    if (class(x) != "reg") {
        stop("x should be a `reg` object.", call. = FALSE)
    }
    result <- x$detail
    return(result)
} 
