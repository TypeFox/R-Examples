#' Subset EEM list
#' 
#' Subset EEM list
#' 
#' @param x EEM class object
#' @param i indices specifying elements to extract
#' @param ... arguments for \code{subset} function
#' 
#' @examples
#' data(applejuice)
#' selected <- applejuice[1-5]
#' 
#' @export
`[.EEM` <- function(x, i, ...) {
    r <- NextMethod("[")
    class(r) <- "EEM"
    r
}