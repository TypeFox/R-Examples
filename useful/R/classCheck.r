#' @title classdf
#' @description Get class information for each column in a \code{\link{data.frame}}.
#' @details Get class information for each column in a \code{\link{data.frame}}.
#' @aliases classdf
#' @export classdf
#' @author Jared P. Lander
#' @param data \code{link{data.frame}} that is to be inspected.
#' @param cols The columns (named or numeric) to be included in the check.
#' @return A vector detailing the class of each column.
#' @examples
#' classdf(CO2)
#' classdf(iris)
#' classdf(mtcars)
classdf <- function(data, cols)
{
    # stop if it is not a data.frame
    if(!is.data.frame(data))
    {
        stop("data must be a data.frame")
    }
    
    # if cols is not supplied check all columns
    if(missing(cols))
    {
        cols <- 1:ncol(data)
    }
    
    # check the class of the columns
    sapply(data[, cols], class)
}