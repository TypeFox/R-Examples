#' @title Penalty of mixture weights
#'
#' @description 
#' Computes the penalty \eqn{\Omega(\mathbf{W})} of the weight matrix (or vector) for a mixture model.
#'
#' @param weigh.matrix \eqn{N \times K} weight matrix
#' @param type type of penalty: \code{c("entropy", "1-Lq", "lognorm")}. 
#' Default: \code{"entropy"}
#' @param q exponent for \eqn{L_q} norm.
#' @param row.average logical; if \code{TRUE} (default) then an average penalty
#' over all rows will be returned (one single number); if \code{FALSE} a vector 
#' of length \eqn{N} will be returned.
#' @param base logarithm base for the \code{"entropy"} penalty. 
#' Default: \code{base = 2}.  Any other real number is allowed; 
#' if \code{base = "num.states"} then it will internally assign it 
#' \code{base = ncol(weigh.matrix)}.
#' @keywords manip array
#' @export
#' @seealso \code{\link{compute_LICORS_loglik}} \code{\link{compute_NEC}}
#' @examples
#' WW = matrix(c(rexp(10, 1/10), runif(10), 1 / 10), ncol=3, byrow = FALSE)
#' WW[1, 1] = 0
#' WW = normalize(WW)
#' compute_mixture_penalty(WW, row.average = FALSE)
#' compute_mixture_penalty(WW, row.average = TRUE) # default: average penalty

compute_mixture_penalty <- function(weigh.matrix, 
                                    type = 
                                      c("entropy", "Lq", "lognorm", "MDL"),
                                    q = 2, row.average = TRUE, base = 2) {
  qq <- q
  num.states <- ncol(weigh.matrix)
  
  type <- match.arg(type)
  if (base == "num.states"){
    base <- num.states
  }

  switch(type, 
         Lq = {
           min.norm <- (1/num.states^(qq - 1))^(1/qq)
           # f(x) = a x + b; b = -a, and a = 1 / (min.norm(w) - 1))
           aa <- 1/(min.norm - 1)
           penalty <- aa * ((rowSums(weigh.matrix^qq))^(1/qq) - 1)
         },
         entropy = {
           temp <- - weigh.matrix * log(weigh.matrix, base)
           penalty <- rowSums(temp, na.rm = TRUE)
         },
         MDL = {
           penalty <- - log(rowSums(weigh.matrix^qq)^(1/qq))
         },
         lognorm = {
           num.samples <- nrow(weigh.matrix)
           penalty <- rep(log(num.states, base), num.samples)
         })
  
  if (row.average) {
    penalty <- mean(penalty)
  }
  return(penalty)
}