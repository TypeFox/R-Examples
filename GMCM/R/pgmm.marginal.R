#' @rdname dgmcm.loglik
pgmm.marginal <- function (z, theta) {
  if (!is.theta(theta)) {
    stop("theta is formatted incorrectly. is.theta(theta) is not TRUE.")
  }
  if (!is.matrix(z)) {
    stop("z is not a matrix. is.matrix(z) is not TRUE.")
  }
  if (ncol(z) != theta$d) {
    stop("Number of columns of z does not equal theta$d.")
  }
  pgmm_marginal(z = z, mus = theta$mu, sigmas = theta$sigma, pie = theta$pie)
}

# pgmm.marginal2 <- function (z, theta) {
#   if (!is.matrix(z))
#     z <- rbind(z)
#   # z is a n by d matrix. (1st col is evaluated in 1st marginal and so on.)
#   tz <- t(z)
#   Margs <- function(k)
#     t(theta$pie[[k]]*pnorm(tz, mean = theta$mu[[k]],
#                            sd = sqrt(diag(theta$sigma[[k]]))))
#   components <- lapply(1:theta$m, Margs)
#   res <- Reduce('+', components)
#
#   return(res)
# }

# pgmm.marginal3 <- function (z, theta) {
#   if (is.matrix(z))
#     z <- cbind(z)
#   # x is a n by d matrix. (1st col is evaluated in 1st marginal and so on.)
#   TempFuncMarginal <- function (k,j) {
#     theta$pie[k]*pnorm(cbind(z)[,j], mean = theta$mu[[k]][j],
#                        sd = sqrt(diag(theta$sigma[[k]])[j]))
#   }
#   TempAggregate <- function (j) {
#     rowSums(rbind(sapply(1:theta$m, FUN = TempFuncMarginal, j = j)))
#   }
#   return(sapply(1:theta$d, FUN = TempAggregate))
# }
