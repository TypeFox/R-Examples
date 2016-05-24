#' Generalized Laguerre polynomials
#' 
#' Generalized Laguerre polynomials as computed by orthopolynom.
#' 
#' @param degree degree of polynomial
#' @param alpha generalization constant
#' @param indeterminate indeterminate
#' @param normalized provide normalized coefficients
#' @return a mpoly object or mpolyList object
#' @author David Kahle calling code from the orthopolynom package
#' @seealso \code{\link{glaguerre.polynomials}}, 
#'   \url{http://en.wikipedia.org/wiki/Laguerre_polynomials}
#' @export
#' @examples
#' 
#' laguerre(0)
#' laguerre(1)
#' laguerre(2)
#' laguerre(3)
#' laguerre(4)
#' laguerre(5)
#' laguerre(6)
#' 
#' laguerre(2)
#' laguerre(2, normalized = TRUE)
#' 
#' laguerre(0:5) 
#' laguerre(0:5, normalized = TRUE)
#' laguerre(0:5, indeterminate = "t")
#' 
#' 
#' 
#' # visualize the laguerre polynomials
#' 
#' library(ggplot2); theme_set(theme_classic())
#' library(tidyr)
#' 
#' s <- seq(-5, 20, length.out = 201)
#' N <- 5 # number of laguerre polynomials to plot
#' (lagPolys <- laguerre(0:N))
#' 
#' # see ?bernstein for a better understanding of
#' # how the code below works
#' 
#' df <- data.frame(s, as.function(lagPolys)(s))
#' names(df) <- c("x", paste0("L_", 0:N))
#' mdf <- gather(df, degree, value, -x)
#' qplot(x, value, data = mdf, geom = "line", color = degree)
#' 
#' qplot(x, value, data = mdf, geom = "line", color = degree) +
#'   coord_cartesian(ylim = c(-25, 25))
#' 
#' 
laguerre <- function(degree, alpha = 0, indeterminate = "x", normalized = FALSE){
  
  stopifnot(all(is.wholenumber(degree)))
  stopifnot(all(degree >= 0))
  
  ## make coefs
  coefs <- glaguerre.polynomials(max(degree), alpha, normalized = normalized)
  
  ## if only one degree is wanted, return that
  if(length(degree) == 1){
    coefs <- rev.default(coefs)[[1]]
    p <- as.mpoly.polynomial(coefs, indeterminate)
    class(p) <- c("laguerre", "mpoly")
    attr(p, "laguerre") <- list(
      degree = length(coefs)-1, 
      alpha = alpha, 
      indeterminate = indeterminate, 
      normalized = normalized
    )
    return(p)
  }
  
  ## if several are wanted, return them
  coefs <- coefs[degree+1]
  ps <- lapply(coefs, function(polynomial){
    p <- as.mpoly.polynomial(polynomial, indeterminate)
    class(p) <- c("laguerre", "mpoly")
    attr(p, "laguerre") <- list(
      degree = length(coefs)-1, 
      alpha = alpha, 
      indeterminate = indeterminate, 
      normalized = normalized
    )
    p
  })
  class(ps) <- "mpolyList"
  ps
  
}



