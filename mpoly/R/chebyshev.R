#' Chebyshev polynomials
#' 
#' Chebyshev polynomials as computed by orthopolynom.
#' 
#' @param degree degree of polynomial
#' @param kind "t" or "u" (Chebyshev polynomials of the first and
#'   second kinds), or "c" or "s"
#' @param indeterminate indeterminate
#' @param normalized provide normalized coefficients
#' @return a mpoly object or mpolyList object
#' @author David Kahle calling code from the orthopolynom package
#' @seealso \code{\link{chebyshev.t.polynomials}}, 
#'   \code{\link{chebyshev.u.polynomials}}, 
#'   \code{\link{chebyshev.c.polynomials}}, 
#'   \code{\link{chebyshev.s.polynomials}}, 
#'   \url{http://en.wikipedia.org/wiki/Chebyshev_polynomials}
#' @export
#' @examples
#' 
#' chebyshev(0)
#' chebyshev(1)
#' chebyshev(2)
#' chebyshev(3)
#' chebyshev(4)
#' chebyshev(5)
#' chebyshev(6)
#' chebyshev(10)
#' 
#' chebyshev(0:5) 
#' chebyshev(0:5, normalized = TRUE)
#' chebyshev(0:5, kind = "u")
#' chebyshev(0:5, kind = "c")
#' chebyshev(0:5, kind = "s")
#' chebyshev(0:5, indeterminate = "t")
#' 
#' 
#' 
#' # visualize the chebyshev polynomials
#' 
#' library(ggplot2); theme_set(theme_classic())
#' library(tidyr)
#' 
#' s <- seq(-1, 1, length.out = 201)
#' N <- 5 # number of chebyshev polynomials to plot
#' (chebPolys <- chebyshev(0:N))
#' 
#' # see ?bernstein for a better understanding of
#' # how the code below works
#' 
#' df <- data.frame(s, as.function(chebPolys)(s))
#' names(df) <- c("x", paste0("T_", 0:N))
#' mdf <- gather(df, degree, value, -x)
#' qplot(x, value, data = mdf, geom = "line", color = degree)
#' 
#' 
#' 
chebyshev <- function(degree, kind = "t", indeterminate = "x", normalized = FALSE){
  
  stopifnot(all(is.wholenumber(degree)))
  stopifnot(all(degree >= 0))

  
  ## deal with kind
  stopifnot(kind %in% c("t","u","c","s"))
  
  ## make coefs
  if(kind == "t") coefs <- chebyshev.t.polynomials(max(degree), normalized)
  if(kind == "u") coefs <- chebyshev.u.polynomials(max(degree), normalized)
  if(kind == "c") coefs <- chebyshev.c.polynomials(max(degree), normalized)
  if(kind == "s") coefs <- chebyshev.s.polynomials(max(degree), normalized)
  
  ## if only one degree is wanted, return that
  if(length(degree) == 1){
    coefs <- rev.default(coefs)[[1]]
    p <- as.mpoly.polynomial(coefs, indeterminate)
    class(p) <- c("chebyshev", "mpoly")
    attr(p, "chebyshev") <- list(
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
    class(p) <- c("chebyshev", "mpoly")
    attr(p, "chebyshev") <- list(
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



