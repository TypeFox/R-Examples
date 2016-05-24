## simplr.R
##
## SimplR - Basic Symbolic Expression Simplification
## 2015 sourcewerk UG
## Oliver Flasch (of@sourcewerk.de)
## Felix Gorschlueter (fg@sourcewerk.de)
## All rights reserved. 
##

##' Basic Symbolic Expression Simplification
##'
##' \code{simplify} simplifies \code{sexp} by applying basic algebraic
##' simplification rules. 
##' \code{simplifyq} quotes its argument, i.e. \code{simplifyq(X)} is
##' equivalent to \code{simplify(quote(X))}.
##'
##' \code{simplify} is a S3 generic method with support for objects of class
##' \code{numeric}, \code{integer}, \code{name}, \code{call}, and
##' \code{function}.
##' SimplR uses code from the Ev3 computer algebra system to implement
##' expression simplification. The following simplification steps are
##' performed:
##' \itemize{
##'   \item consolidate product coefficients 
##'   \item distribute coefficients over sums  
##'   \item convert differences to sums
##'   \item simplify constants
##'   \item simplify products
##'   \item compact linear parts
##'   \item simplify trigonometrics
##' }
##'
##' @param sexp An R object to simplify. See details.
##' @return The simplified expression. 
##'
##' @examples
##' simplifyq(3*2+1)              #=> 7
##'
##' simplifyq(1 * x)              #=> x
##' simplifyq(x / x)              #=> 1
##' simplifyq(x - x)              #=> 0
##' simplifyq(x + 1 - 1)          #=> x
##'
##' simplifyq(f(x) + f(x) + y)    #=> y + 2 * f(x)
##' simplifyq(sin(x)^2+cos(x)^2)  #=> 1
##'
##' simplify(function(a,b) a + a + 3 * f(b) * 5 / f(b))
##' #=> function(a, b) 15 + 2 * a
##'
##' @useDynLib simplr
##' @rdname simplify
##' @seealso \url{http://www.lix.polytechnique.fr/~liberti/Ev3.pdf}
##' @export
simplify <- function(sexp) UseMethod("simplify")

##' @rdname simplify
##' @export
simplifyq <- function(sexp) simplify(substitute(sexp))
 
##' @rdname simplify
##' @method simplify call 
##' @export
simplify.call <- function(sexp) {
  .Call("simplify", sexp, PACKAGE = "simplr")
}

##' @rdname simplify
##' @method simplify function 
##' @export
simplify.function <- function(sexp) {
  simplifiedFun <- sexp # copy
  body(simplifiedFun) <- simplify.call(body(sexp))
  return (simplifiedFun)
}

##' @rdname simplify
##' @method simplify numeric 
##' @export
simplify.numeric <- simplify.call

##' @rdname simplify
##' @method simplify integer 
##' @export
simplify.integer <- simplify.call

##' @rdname simplify
##' @method simplify name 
##' @export
simplify.name <- simplify.call

