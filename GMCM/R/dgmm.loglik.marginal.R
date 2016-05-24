#' @rdname dgmcm.loglik
dgmm.loglik.marginal <- function (theta, x, marginal.loglik = TRUE) {
  if (!is.theta(theta))
    stop("theta is not formatted correctly")
  if (!is.matrix(x))
    stop("x is not a matrix")
  if (ncol(x) != theta$d)
    stop("Number of colums of x does not equal theta$d")

  dgmm_loglik_marginal(mus = theta$mu,
                       sigmas = theta$sigma,
                       pie = theta$pie,
                       z = x,
                       marginal_loglik = marginal.loglik)
}

# dgmm.loglik.marginal2 <- function (theta, x, marginal.loglik = TRUE) {
#   TempFuncMarginal <- function (k, j) {
#     theta$pie[k]*dnorm(cbind(x)[, j], mean = theta$mu[[k]][j],
#                        sd = sqrt(diag(theta$sigma[[k]])[j]))
#   }
#   TempAggregate <- function (j) {
#     rowSums(rbind(sapply(1:theta$m, FUN = TempFuncMarginal, j = j)))
#   }
#   loglik <- log(sapply(1:ncol(cbind(x)), FUN = TempAggregate))
#   if (!marginal.loglik)
#     loglik <- colSums(loglik)
#   return(loglik)
# }
#
# dgmm.loglik.marginal3 <- function(theta, x, marginal.loglik = TRUE) {
#   loglik <- NULL
#   for (j in 1:ncol(cbind(x))) {
#     tmp <- 0
#     for (k in 1:theta$m ) {
#       tmp <- tmp +
#         theta$pie[k]*dnorm(cbind(x)[,j],
#                            mean = theta$mu[[k]][j],
#                            sd = sqrt(diag(theta$sigma[[k]])[j]))
#       #print(cbind(theta$pie[k]*dnorm(cbind(x)[,j],
#       #                         mean = theta$mu[[k]][j],
#       #                         sd = sqrt(diag(theta$sigma[[k]])[j]))))
#     }
#     loglik <- cbind(loglik, log(tmp))
#   }
#
#
#   if (!marginal.loglik)
#     loglik <- colSums(loglik)
#
#   return(loglik)
# }
