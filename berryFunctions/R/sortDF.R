#' sort dataframes by column
#' 
#' sort a data.frame by column - basically just a wrapper for order
#' 
#' @return data.frame
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, June 2015
#' @seealso \code{\link{sort}}, \code{\link{order}}
#' @keywords univar manip arith
#' @export
#' @examples
#' 
#' sortDF(USArrests[USArrests$Murder>8,], "Assault")
#' sortDF(USArrests[USArrests$Murder>8,], 3, decreasing=TRUE)
#' 
#' @param df Data.frame to be sorted
#' @param col Column (index or name) to be sorted by
#' @param \dots Further arguments passed to \code{\link{order}}, like eg \code{decreasing=TRUE}
#' 
sortDF <- function(
df,
col,
...
)
{
df[order(df[,col], ...),]
}
