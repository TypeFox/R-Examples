#' Find the minimum and maximum of a vector
#' 
#' The function \code{mm()} finds the minimum and maximum of a vector. It is 
#' intended for use with \code{eplot()} to properly scale the axes.
#' 
#' 
#' @aliases mm
#' @param x a vector
#@note %% ~~further notes~~
#' @author Carlisle Rainey (\href{mailto:carlislerainey@@gmail.com}{e-mail},
#' \href{http://www.carlislerainey.com}{website})
#@seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#@references %% ~put references to the literature/web site here ~
#@keywords ~kwd1 ~kwd2
#' @examples
#' 
#' x <- rnorm(100)
#' y <- rnorm(100)
#' 
#' par(mfrow = c(1,1), mar = c(5,4,4,2), oma = c(0,0,0,0))
#' eplot(x, y, xlim = mm(x), ylim = mm(y))
#' 
#' @export mm
mm <- function(x) {
  min <- min(x, na.rm = TRUE)
  max <- max(x, na.rm = TRUE)
  res <- c(min, max)
  return(res)
}