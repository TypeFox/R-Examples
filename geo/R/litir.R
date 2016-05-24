#' Display palette in effect
#' 
#' A grid of colors is plotted.
#' 
#' Palette in effect gets repeated for \code{n} greater than
#' \code{length(palette())}.
#' 
#' @param n Number of columns and rows of colors to display
#' @return No value returned, plots an \code{n} by \code{n} grid of colors.
#' @note Simpler version of \code{\link{Rlitir}}, one of those should do,
#' possibly with a name change.
#' @seealso \code{\link{Rlitir}}
#' @keywords colors
#' @export litir
litir <-
function(n)
{
        x <- c(1:(n + 1))
         plot(x, x)
        y <- x
        for(j in 1:(n - 1)) {
                for(i in 1:(n + 1)) {
                        polygon(c(x[i], x[i + 1], x[i + 1], x[i]), c(y[j],
                                y[j], y[j + 1], y[j + 1]), col = ((j - 1) *
                                n + i - 1))
                        lines(c(x[i], x[i + 1], x[i + 1], x[i], x[i]), c(y[
                                j], y[j], y[j + 1], y[j + 1], y[j]), col = 1)
                        text((x[i] + x[i + 1])/2, (y[j] + y[j + 1])/2,
                                as.character((j - 1) * n + i - 1))
                }
        }
        return(invisible())
}

