#' Convenience function for Gauss-Hermite quadrature
#' 
#' Convenience function for evaluation of Gauss-Hermite quadrature
#' 
#' This function performs classical unidimensional Gauss-Hermite quadrature
#' with the function f using the rule provided; that is, it approximates
#' \deqn{\int_{-\infty}^{\infty} f(x) \exp(-x^2) \, dx}{ integral( f(x)
#' exp(-x^2), -Inf, Inf)} by evaluating \deqn{ \sum_i w_i f(x_i) }{sum( w *
#' f(x) )}
#' 
#' @param f Function to integrate with respect to first (scalar) argument; this
#' does not include the weight function \code{exp(-x^2)}
#' @param rule Gauss-Hermite quadrature rule to use, as produced by
#' \code{\link{gaussHermiteData}}
#' @param ... Additional arguments for f
#' @return Numeric (scalar) with approximation integral of f(x)*exp(-x^2) from
#' -Inf to Inf.
#' @author Alexander W Blocker \email{ablocker@@gmail.com}
#' @export
#' @seealso \code{\link{gaussHermiteData}}, \code{\link{ghQuad}}
#' @references Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss
#' Quadrature Rules. Mathematics of Computation 23 (106): 221-230.
#' 
#' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite Quadrature.
#' Biometrika, 81(3) 624-629.
#' @keywords math
#' @examples
#' 
#' # Get quadrature rules
#' rule10  <- gaussHermiteData(10)
#' rule100 <- gaussHermiteData(100)
#' 
#' # Check that rule is implemented correctly
#' f <- function(x) rep(1,length(x))
#' if (!isTRUE(all.equal(sqrt(pi), ghQuad(f, rule10), ghQuad(f, rule100)))) {
#'   print(ghQuad(f, rule10))
#'   print(ghQuad(f, rule100))
#' }
#' # These should be 1.772454
#' 
#' f <- function(x) x
#' if (!isTRUE(all.equal(0.0, ghQuad(f, rule10), ghQuad(f, rule100)))) {
#'   print(ghQuad(f, rule10))
#'   print(ghQuad(f, rule100))
#' }
#' # These should be zero
#' 
#' 
ghQuad <- function(f, rule, ...) {
    # Integrate function according to given quadrature rule
    # Simple wrapper
    sum(rule$w * f(rule$x, ...))
}



#' Adaptive Gauss-Hermite quadrature using Laplace approximation
#' 
#' Convenience function for integration of a scalar function g based upon its
#' Laplace approximation.
#' 
#' This function approximates \deqn{\int_{-\infty}^{\infty} g(x) \, dx}{
#' integral( g(x), -Inf, Inf)} using the method of Liu & Pierce (1994). This
#' technique uses a Gaussian approximation of g (or the distribution component
#' of g, if an expectation is desired) to "focus" quadrature around the
#' high-density region of the distribution. Formally, it evaluates: \deqn{
#' \sqrt{2} \hat{\sigma} \sum_i w_i \exp(x_i^2) g(\hat{\mu} + \sqrt{2} }{
#' sqrt(2) * sigmaHat * sum( w * exp(x^2) * g(muHat + sqrt(2) * sigmaHat * x))
#' }\deqn{\hat{\sigma} x_i) }{ sqrt(2) * sigmaHat * sum( w * exp(x^2) * g(muHat
#' + sqrt(2) * sigmaHat * x)) } where x and w come from the given rule.
#' 
#' This method can, in many cases (where the Gaussian approximation is
#' reasonably good), achieve better results with 10-100 quadrature points than
#' with 1e6 or more draws for Monte Carlo integration. It is particularly
#' useful for obtaining marginal likelihoods (or posteriors) in hierarchical
#' and multilevel models --- where conditional independence allows for
#' unidimensional integration, adaptive Gauss-Hermite quadrature is often
#' extremely effective.
#' 
#' @param g Function to integrate with respect to first (scalar) argument
#' @param muHat Mode for Laplace approximation
#' @param sigmaHat Scale for Laplace approximation (\code{sqrt(-1/H)}, where H
#' is the second derivative of g at muHat)
#' @param rule Gauss-Hermite quadrature rule to use, as produced by
#' \code{\link{gaussHermiteData}}
#' @param ... Additional arguments for g
#' @return Numeric (scalar) with approximation integral of g from -Inf to Inf.
#' @author Alexander W Blocker \email{ablocker@@gmail.com}
#' @export
#' @seealso \code{\link{gaussHermiteData}}, \code{\link{ghQuad}}
#' @references Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite
#' Quadrature. Biometrika, 81(3) 624-629.
#' @keywords math
#' @examples
#' 
#' # Get quadrature rules
#' rule10  <- gaussHermiteData(10)
#' rule100 <- gaussHermiteData(100)
#' 
#' # Estimating normalizing constants
#' g <- function(x) 1/(1+x^2/10)^(11/2) # t distribution with 10 df
#' aghQuad(g, 0, 1.1,  rule10)
#' aghQuad(g, 0, 1.1,  rule100)
#' # actual is
#' 1/dt(0,10)
#' 
#' # Can work well even when the approximation is not exact
#' g <- function(x) exp(-abs(x)) # Laplace distribution
#' aghQuad(g, 0, 2,  rule10)
#' aghQuad(g, 0, 2,  rule100)
#' # actual is 2
#' 
#' # Estimating expectations
#' # Variances for the previous two distributions
#' g <- function(x) x^2*dt(x,10) # t distribution with 10 df
#' aghQuad(g, 0, 1.1,  rule10)
#' aghQuad(g, 0, 1.1,  rule100)
#' # actual is 1.25
#' 
#' # Can work well even when the approximation is not exact
#' g <- function(x) x^2*exp(-abs(x))/2 # Laplace distribution
#' aghQuad(g, 0, 2,  rule10)
#' aghQuad(g, 0, 2,  rule100)
#' # actual is 2
#' 
#' 
aghQuad <- function(g, muHat, sigmaHat, rule, ...) {
    # Adaptive Gauss-Hermite quadrature as in Liu & Pierce (1994)
    
    # Get transformed nodes
    z <- muHat + sqrt(2)*sigmaHat*rule$x

    # Transform weights to account for use of importance-sampling type ratio
    wStar <- exp(rule$x*rule$x + log(rule$w))
    wStar
    
    # Approximate integrate
    val <- sqrt(2)*sigmaHat*sum(wStar*g(z, ...))

    return(val)
}
