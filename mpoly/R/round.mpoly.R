#' Round the coefficients of a polynomial
#'
#' Round the coefficients of an mpoly object.
#' 
#' @param x an mpoly object
#' @param digits number of digits to round to
#' @return the rounded mpoly object
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @seealso \code{\link{mp}}
#' @export
#' @examples
#' 
#' p <- mp("x + 3.14159265")^4
#' p
#' round(p)
#' round(p, 0)
#' 
#' \dontrun{
#' library(plyr)
#' library(ggplot2)
#' library(stringr)
#' 
#' n <- 101
#' s <- seq(-5, 5, length.out = n)
#' 
#' # one dimensional case
#' df <- data.frame(x = s)
#' df <- mutate(df, y = -x^2 + 2*x - 3 + rnorm(n, 0, 2))
#' qplot(x, y, data = df)
#' mod <- lm(y ~ x + I(x^2), data = df)
#' p <- as.mpoly(mod)
#' qplot(x, y, data = df) +
#'   stat_function(fun = as.function(p), colour = 'red')
#' 
#' p
#' round(p, 1)
#' qplot(x, y, data = df) +
#'   stat_function(fun = as.function(p), colour = 'red') +
#'   stat_function(fun = as.function(round(p,1)), colour = 'blue')
#'
#'
#' }
#' 
round.mpoly <- function(x, digits = 3){
  
  ## round coefficients
  p <- lapply(x, function(term){
    term["coef"] <- round(term["coef"], digits = digits)
    term
  })
  
  ## drop zero terms
  p <- Filter(function(v) v[["coef"]] != 0, p)
  
  ## class
  class(p) <- "mpoly"
  
  ## out
  p
}

