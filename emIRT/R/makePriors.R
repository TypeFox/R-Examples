#' Create an R list representing prior beliefs over parameters.
#'
#' Create an R list representing diffuse prior beliefs over the different
#' parameters in the model. There is little customization possible through this
#' function.
#'
#' @param .D (integer) number of spatial dimensions in model
#' @param .J (integer) number of "bills" in model
#' @param .N (integer) number of "legislators" in model
#'
#' @return A list with elements 'x' and 'beta'. 'x' (or 'beta') contains two
#' elements: 'mu' and 'sigma'. 'mu' is a '.D' by 1 (or '.D+1' by 1) matrix
#' representing the means of the prior distributions. Each of the '.N' 'x'
#' vectors (or each of the 'J' 'beta' vectors) is given the same prior
#' mean. 'sigma' is a 'D' by 'D' (or 'D+1' by 'D+1') matrix representing the
#' covariance matrix of the prior distribution. Each of the '.N' 'x' vectors (or
#' each of the '.J' 'beta' vectors) is given the same prior covariance matrix.
#'
#' @seealso \code{\link{fastest}} which uses the list of priors

makePriors <- function(.N = 20,
                       .J = 100,
                       .D = 1
                       ) {
    out <- vector(mode = "list")

    out$x <- list(mu = matrix(rep(0, .D), nrow = .D),
                  sigma = 1^2 * diag(.D)
                  )

    out$beta <- list(mu = matrix(rep(0, .D + 1), nrow = .D + 1),
                     sigma = 5^2 * diag(.D + 1)
                     )
    return(out)
}


makePriorsE <- function(.N = 20,
                        .J = 100
                        ) {
    out <- vector(mode = "list")

    out$theta <- list(mu = matrix(rep(0, 1), nrow = 1),
                      sigma = 1^2 * diag(1)
                      )

    out$w <- list(mu = matrix(rep(0, 1), nrow = 1),
                  sigma = 1^2 * diag(1)
                  )

    out$alpha <- list(mu = matrix(rep(0, 1), nrow = 1),
                      sigma = 5^2 * diag(1)
                      )

    out$beta <- list(mu = matrix(rep(0, 1), nrow = 1),
                     sigma = 5^2 * diag(1)
                     )

    out$gamma <- list(mu = matrix(rep(1, 1), nrow = 1),
                      sigma = 1^2 * diag(1)
                      )

    return(out)
}
