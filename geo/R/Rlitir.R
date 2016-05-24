#' Display colors.
#' 
#' A grid of colors is plotted.
#' 
#' 
#' @param n Number of columns and rows of colors to display
#' @param col A vector of colors
#' @return No value returned, plots an \code{n} by \code{n} grid of colors.
#' @seealso \code{\link{colorRampPalette}}
#' @keywords color
#' @examples
#' 
#' # simple, perhaps not so useful application with default palette:
#' 
#' Rlitir(12, 1:144)
#' 
#' # Define a palette with some colors:
#' 
#' ramp <- colorRampPalette(c("khaki1", "gold", "orange", 
#'   "darkorange2", "red", "darkred", "black"))
#' 
#' # number of columns and rows to display
#' 
#' n <- 10
#' 
#' Rlitir(n, ramp(n^2))
#' 
#' @export Rlitir
Rlitir <-
function(n,col)
{
        x <- c(1:(n + 1))
         plot(x, x)
        y <- x
        for(j in 1:(n - 1)) {
                for(i in 1:(n + 1)) {
                        polygon(c(x[i], x[i + 1], x[i + 1], x[i]), c(y[j],
                                y[j], y[j + 1], y[j + 1]), col = col[((j - 1) *
                                n + i - 1)])
                        lines(c(x[i], x[i + 1], x[i + 1], x[i], x[i]), c(y[
                                j], y[j], y[j + 1], y[j + 1], y[j]), col = 1)
                        text((x[i] + x[i + 1])/2, (y[j] + y[j + 1])/2,
                                as.character((j - 1) * n + i - 1))
                }
        }
        return(invisible())
}

