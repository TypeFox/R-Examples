#' @rdname beta-utils
#' @description
#' \code{check_beta} checks if \eqn{\boldsymbol \beta} defines a
#'     valid distribution, e.g., for normal distribution \code{'sigma'} must be
#'     positive.
#' @return \code{check_beta} throws an error if \eqn{\boldsymbol \beta} is not
#'     appropriate for the given distribution; e.g., if it has too many values
#'     or if they are not within proper bounds (e.g., \code{beta['sigma']} of a
#'     \code{"normal"} distribution must be positive).
#' @export
#' @examples
#' 
#' \dontrun{
#' check_beta(beta = c(1, 1, -1), distname = "normal")
#' }
#' 
check_beta <- function(beta, distname) {

  stopifnot(is.numeric(beta),
            !anyNA(beta))
  check_distname(distname)
  dist.text <- paste0("For a '", distname, "' distribution 'beta' must ")
  switch(distname, 
         "normal" = {
           if (length(beta) != 2) {
             stop(dist.text, " be a vector of length 2.",
                  "\n Currently it has length ", length(beta), ".")
           } 
           if (beta[2] <= 0) {
             stop(dist.text, " must have a non-negative second entry (beta[2] ~ standard deviation).", 
                  "Currently, beta[2] = ", beta[2])
           }
         },
         "t" = {
           if (length(beta) != 3) {
             stop(dist.text, " be a vector of length 3.",
                  "\n Currently it has length ", length(beta), ".")
           } 
           if (beta[2] < 0) {
             stop(dist.text, " must have a positive second entry (beta[2] ~ scale parameter).", 
                  "Currently, beta[2] = ", beta[2])
           }
           if (beta[3] <= 0) {
             stop(dist.text, " must have a positive third entry (beta[3] ~ degrees of freedom).",
                  "Currently, beta[3] = ", beta[3])
           }
         },
        "exp" = {
          if (length(beta) != 1) {
             stop(dist.text, " must be a scalar or vector of length 1.") 
           }
           if (beta <= 0) {
             stop(dist.text, " must be positive.  It is the rate of the exponential distribution.",
                  "Currently, beta = ", beta) 
           }
         },
        "f" = {
          if (length(beta) != 2) {
            stop(dist.text, " be a vector of length 2 (lower and upper bound).",
                 "\n Currently it has length ", length(beta), ".")
          } 
          if (any(beta < 0)) {
            stop(dist.text, " must be positive.  It is the rate of the exponential distribution.",
                 "Currently, beta[1] = ", beta[1], " and beta[2] = ", beta[2])
          } 
          if (beta[2] <= 4) {
            stop(dist.text, " must have finite variance; thus df2 > 4.",
                 "Currently, df2 = ", beta[2], ".")
          }
        },
         "unif" = {
           if (length(beta) != 2) {
             stop(dist.text, " be a vector of length 2 (lower and upper bound).",
                  "\n Currently it has length ", length(beta), ".")
           } 
           if (beta[1] > beta[2]) {
             stop(dist.text, " must have beta[1] < beta[2], since beta[1] = lower bound and beta[2] = upper bound",
                  " of the uniform distribution U(a, b).", 
                  "Currently, beta[1] = ", beta[1], " and beta[2] = ", beta[2])
           }
         })
}
