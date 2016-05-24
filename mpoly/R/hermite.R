#' Hermite polynomials
#' 
#' Hermite polynomials as computed by orthopolynom.
#' 
#' @param degree degree of polynomial
#' @param kind "he" (default, probabilists', see Wikipedia article)
#'   or "h" (physicists)
#' @param indeterminate indeterminate
#' @param normalized provide normalized coefficients
#' @return a mpoly object or mpolyList object
#' @author David Kahle calling code from the orthopolynom package
#' @seealso \code{\link{hermite.h.polynomials}}, 
#'   \code{\link{hermite.he.polynomials}}, 
#'   \url{http://en.wikipedia.org/wiki/Hermite_polynomials}
#' @export
#' @examples
#' 
#' hermite(0)
#' hermite(1)
#' hermite(2)
#' hermite(3)
#' hermite(4)
#' hermite(5)
#' hermite(6)
#' hermite(10)
#' 
#' hermite(0:5) 
#' hermite(0:5, normalized = TRUE)
#' hermite(0:5, indeterminate = "t")
#' 
#' 
#' 
#' # visualize the hermite polynomials
#' 
#' library(ggplot2); theme_set(theme_classic())
#' library(tidyr)
#' 
#' s <- seq(-3, 3, length.out = 201)
#' N <- 5 # number of hermite polynomials to plot
#' (hermPolys <- hermite(0:N))
#' 
#' # see ?bernstein for a better understanding of
#' # how the code below works
#' 
#' df <- data.frame(s, as.function(hermPolys)(s))
#' names(df) <- c("x", paste0("T_", 0:N))
#' mdf <- gather(df, degree, value, -x)
#' qplot(x, value, data = mdf, geom = "line", color = degree)
#' 
#' 
#' 
hermite <- function(degree, kind = "he", indeterminate = "x", normalized = FALSE){
  
  stopifnot(all(is.wholenumber(degree)))
  stopifnot(all(degree >= 0))
  
  
  ## deal with kind
  stopifnot(kind %in% c("h","he"))
  
  ## make coefs
  if(kind == "h")  coefs <-  hermite.h.polynomials(max(degree), normalized)
  if(kind == "he") coefs <- hermite.he.polynomials(max(degree), normalized)

  
  ## if only one degree is wanted, return that
  if(length(degree) == 1){
    coefs <- rev.default(coefs)[[1]]
    p <- as.mpoly.polynomial(coefs, indeterminate)
    class(p) <- c("hermite", "mpoly")
    attr(p, "hermite") <- list(
      degree = length(polynomial)-1, 
      kind = kind, 
      indeterminate = indeterminate,
      normalized = normalized
    )
    return(p)
  }
  
  ## if several are wanted, return them
  coefs <- coefs[degree+1]
  ps <- lapply(coefs, function(polynomial){
    p <- as.mpoly.polynomial(polynomial, indeterminate)
    class(p) <- c("hermite", "mpoly")
    attr(p, "hermite") <- list(
      degree = length(polynomial)-1, 
      kind = kind, 
      indeterminate = indeterminate,
      normalized = normalized
    )
    p
  })
  class(ps) <- "mpolyList"
  ps
  
}



