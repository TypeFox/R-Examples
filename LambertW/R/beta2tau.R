#' @rdname beta-utils
#' @description
#' 
#' \code{beta2tau} converts \eqn{\boldsymbol \beta} to the transformation vector
#'     \eqn{\tau = (\mu_x, \sigma_x, \gamma = 0, \alpha = 1, \delta = 0)}, which
#'     defines the Lambert W\eqn{\times} F random variable mapping from \eqn{X}
#'     to \eqn{Y} (see \code{\link{tau-utils}}). Parameters \eqn{\mu_x} and
#'     \eqn{\sigma_x} of \eqn{X} in general depend on \eqn{\boldsymbol \beta}
#'     (and may not even exist for \code{use.mean.variance = TRUE}; in this case
#'     \code{beta2tau} will throw an error).
#' 
#' @return
#' \code{beta2tau} returns a numeric vector, which is \eqn{\tau =
#'     \tau(\boldsymbol \beta)} implied by \code{beta} and \code{distname}.
#' @export
#' @seealso
#' \code{\link{tau-utils}}, \code{\link{theta-utils}}
#' @examples
#' # By default: delta = gamma = 0 and alpha = 1
#' beta2tau(c(1, 1), distname = "normal") 
#' \dontrun{
#'   beta2tau(c(1, 4, 1), distname = "t")
#' }
#' beta2tau(c(1, 4, 1), distname = "t", use.mean.variance = FALSE)
#' beta2tau(c(1, 4, 3), distname = "t") # no problem
#' 
beta2tau <- function(beta, distname, use.mean.variance = TRUE) {
  stopifnot(is.numeric(beta),
            is.logical(use.mean.variance))
  
  check_distname(distname)
  check_beta(beta, distname = distname)
  names(beta) <- get_beta_names(distname)
  # sigma_x is standard deviation, not variance.
  tau <- c()

  switch(distname,
         chisq = {
           if (use.mean.variance) {
             tau[c("mu_x", "sigma_x")] <- c(0, 2 * beta)
           } else {
             tau[c("mu_x", "sigma_x")] <- c(0, beta)
           }
         },
         exp = {
           # scale = std dev
           tau[c("mu_x", "sigma_x")] <- c(0, 1 / beta[1])
         },
         "f" = {
           if (use.mean.variance) {
             tau["mu_x"] <- 0
             tau["sigma_x"] <- sqrt((2 * beta[2]^2 * (beta[1] + beta[2] - 2))/
                                      (beta[1] * (beta[2] - 2)^2 * (beta[2] - 4)))
           } else {
             tau[c("mu_x", "sigma_x")] <- c(0, 1)
           }
          },
         gamma = {
           if (use.mean.variance) {
             tau["mu_x"] <- 0
             tau["sigma_x"] <- sqrt(beta["shape"]) * beta["scale"]
           } else {
             tau["mu_x"] <- 0
             tau["sigma_x"] <- beta["scale"]
           }
         },
         laplace = {
             tau["mu_x"] <- beta[1]
             tau["sigma_x"] <- sqrt(2) * beta[2]
         },
         normal = {
           # mean = location and scale = std dev
           tau[c("mu_x", "sigma_x")] <- beta[1:2]
         },
         t = {
           if (use.mean.variance) {
             nu <- beta[3]
             if (nu <= 2) {
               stop("A t-distribution with df = ", nu,
                    " does not have finite variance, which is required",
                    " for a mean-variance Lambert W x t distribution.",
                    " If your data seems to have non-finite second moments ",
                    "consider using a general location-scale Lambert W x F ",
                    " distribution by setting 'use.mean.variance = FALSE'.")
             }
             ss <- beta[2]
             scaling.factor <- sqrt(nu/(nu - 2))
             tau["mu_x"] <- beta[1]
             tau["sigma_x"] <- ss * scaling.factor
           } else {
             tau[c("mu_x", "sigma_x")] <- beta[c("location", "scale")] 
           }
         },
         unif = {  
           if (use.mean.variance) {
             tau["mu_x"] <- 0.5 * (beta[1] + beta[2])
             tau["sigma_x"] <- sqrt(1/12 * (beta[2] - beta[1])^2)
           } else {
             if (beta[1] != 0) {
                warning("The lower limit of the uniform distribution is ",
                        " non-zero.  When using the general location-scale",
                        " version of the Lambert W x F distribution",
                        " the uniform distribution is usually of scale type.")
             }
             tau["mu_x"] <- 0
             tau["sigma_x"] <- beta[2]
           }
         }
  )
  if (length(tau) == 0) {
    stop("Seems like distribution '", distname, "' is not supported.")    
  }

  # use default values here
  tau <- c(tau, c(alpha = 1, gamma = 0, delta = 0))
  check_tau(tau)
  return(tau)
} 
