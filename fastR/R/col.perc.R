#' Row and Column Percentages
#' 
#' Convenience wrappers around \code{\link{apply}()} to compute row and column
#' percentages of matrix-like structures, including output of
#' \code{\link{xtabs}}.
#' 
#' 
#' @aliases col.perc row.perc
#' @param x matrix-like structure
#' @author Randall Pruim
#' @keywords array
#' @export
#' @examples
#' 
#' row.perc(xtabs(~Airline+Result,airlineArrival))
#' 
col.perc <- function (x) 
{
    x/rep(apply(x, 2, sum), each = dim(x)[1])
}

#' @rdname col.perc
#' @export
row.perc <- function (x) 
{
    x/apply(x, 1, sum)
}
