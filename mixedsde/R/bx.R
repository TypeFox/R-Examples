#' Computation Of The Drift Coefficient

#' @description Computation of the drift coefficient
#' @param x vector of data
#' @param fixed drift constant in front of X (when there is one additive random effect), 0 otherwise 
#' @param random 1 if there is one additive random effect, 2 one multiplicative random effect or c(1,2) for 2 random effects 
#' @return
#' \item{b}{The drift is \eqn{b(x,\phi) = \phi_1 b_1(x) + \phi_2 b_2(x)}, the output is \eqn{b_2} except when random c(1,2) then the output is the vector \eqn{(b_1,b_2)^t}}
#' @keywords drift



bx <- function(x, fixed, random) {
    
    if (sum(random) > 2) {
        return(matrix(c(1, -x), 2, 1))
    }
    if (sum(random) == 1) {
        return(-fixed * x)
    }
    if (sum(random) == 2) {
        return(-x)
    }
} 
