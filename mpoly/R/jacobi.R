#' Jacobi polynomials
#' 
#' Jacobi polynomials as computed by orthopolynom.
#' 
#' @param degree degree of polynomial
#' @param alpha the first parameter, also called p
#' @param beta the second parameter, also called q
#' @param kind "g" or "p"
#' @param indeterminate indeterminate
#' @param normalized provide normalized coefficients
#' @return a mpoly object or mpolyList object
#' @author David Kahle calling code from the orthopolynom package
#' @seealso \code{\link{jacobi.g.polynomials}}, 
#'   \code{\link{jacobi.p.polynomials}} 
#'   \url{http://en.wikipedia.org/wiki/Jacobi_polynomials}
#' @export
#' @examples
#' 
#' jacobi(0)
#' jacobi(1)
#' jacobi(2)
#' jacobi(3)
#' jacobi(4)
#' jacobi(5)
#' jacobi(6)
#' jacobi(10, 2, 2, normalized = TRUE)
#' 
#' jacobi(0:5) 
#' jacobi(0:5, normalized = TRUE)
#' jacobi(0:5, kind = "g")
#' jacobi(0:5, indeterminate = "t")
#' 
#' 
#' 
#' # visualize the jacobi polynomials
#' 
#' library(ggplot2); theme_set(theme_classic())
#' library(tidyr)
#' 
#' s <- seq(-1, 1, length.out = 201)
#' N <- 5 # number of jacobi polynomials to plot
#' (jacPolys <- jacobi(0:N, 2, 2))
#' 
#' df <- data.frame(s, as.function(jacPolys)(s))
#' names(df) <- c("x", paste0("P_", 0:N))
#' mdf <- gather(df, degree, value, -x)
#' qplot(x, value, data = mdf, geom = "line", color = degree)
#'   
#' qplot(x, value, data = mdf, geom = "line", color = degree) +
#'   coord_cartesian(ylim = c(-30, 30))
#' 
#' 
#' 
jacobi <- function(degree, alpha = 1, beta = 1, kind = "p", indeterminate = "x", normalized = FALSE){
  
  stopifnot(all(is.wholenumber(degree)))
  stopifnot(all(degree >= 0))
  
  
  ## deal with kind
  stopifnot(kind %in% c("g","p"))
  
  ## make coefs
  if(kind == "g") coefs <- jacobi.g.polynomials(max(degree), p = alpha, q = beta, normalized)
  if(kind == "p") coefs <- jacobi.p.polynomials(max(degree), alpha = alpha, beta = beta, normalized)
  
  ## if only one degree is wanted, return that
  if(length(degree) == 1){
    coefs <- rev.default(coefs)[[1]]
    p <- as.mpoly.polynomial(coefs, indeterminate)
    class(p) <- c("jacobi", "mpoly")
    attr(p, "jacobi") <- list(
      degree = length(polynomial)-1, 
      kind = kind, 
      indeterminate = indeterminate,
      normalized = normalized,
      alpha = alpha, beta = beta
    )
    return(p)
  }
  
  ## if several are wanted, return them
  coefs <- coefs[degree+1]
  ps <- lapply(coefs, function(polynomial){
    p <- as.mpoly.polynomial(polynomial, indeterminate)
    class(p) <- c("jacobi", "mpoly")
    attr(p, "jacobi") <- list(
      degree = length(polynomial)-1, 
      kind = kind, 
      indeterminate = indeterminate,
      normalized = normalized,
      alpha = alpha, beta = beta
    )
    p
  })
  class(ps) <- "mpolyList"
  ps
  
}



