#' Goodman and Kruskal's lambda
#'
#' Computes Goodman and Kruskal's lambda for given table.
#' @param table a \code{table} of two variables or a \code{data.frame} representation of the cross-table of the two variables without marginals
#' @param direction numeric value of \code{c(0, 1, 2)} where 1 means the lambda value computed for row, 2 for columns and 0 for both
#' @return computed lambda value(s) for row/col of given table
#' @examples \dontrun{
#' ## quick example
#' x <- data.frame(x = c(5, 4, 3), y = c(9, 8, 7), z = c(7, 11, 22), zz = c(1, 15, 8))
#' lambda.test(x)    # 0.1 and 0.18333
#' lambda.test(t(x)) # 0.18333 and 0.1
#'
#' ## historical data (see the references above: p. 744)
#' men.hair.color <- data.frame(
#' b1 = c(1768, 946, 115),
#' b2 = c(807, 1387, 438),
#' b3 = c(189, 746, 288),
#' b4 = c(47, 53, 16)
#' )
#' row.names(men.hair.color) <- paste0('a', 1:3)
#' lambda.test(men.hair.color)
#' lambda.test(t(men.hair.color))
#'
#' ## some examples on mtcars
#' lambda.test(table(mtcars$am, mtcars$gear))
#' lambda.test(table(mtcars$gear, mtcars$am))
#' lambda.test(table(mtcars$am, mtcars$gear), 1)
#' lambda.test(table(mtcars$am, mtcars$gear), 2)
#' }
#' @references \itemize{
#'   \item Goodman, L.A., Kruskal, W.H. (1954) Measures of association for cross classifications. Part I. \emph{Journal of the American Statistical Association} \bold{49}, 732--764
#' }
#' @export
lambda.test <- function(table, direction = 0) {

    if (!is.numeric(direction))
        stop('Direction should be an integer between 0 and 2!')
    if (!direction %in% c(0, 1, 2))
        stop('Direction should be an integer between 0 and 2!')

    if (direction != 0)
        (base::sum(as.numeric(apply(table, direction, base::max))) - ifelse(direction == 1, base::max(colSums(table)), base::max(rowSums(table)))) / (base::sum(table) - ifelse(direction == 1, base::max(colSums(table)), base::max(rowSums(table))))
    else
        list(row = lambda.test(table, 1), col = lambda.test(table, 2))

}
