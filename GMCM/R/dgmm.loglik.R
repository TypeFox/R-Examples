#' @rdname dgmcm.loglik
dgmm.loglik <- function (theta, z, marginal.loglik = FALSE) {  
  if (!is.theta(theta))
    stop("theta is not formatted correctly")
  if (!is.matrix(z))
    stop("z is not a matrix")
  if (ncol(z) != theta$d)
    stop("Number of colums of z does not equal theta$d")
  
  dgmm_loglik(mus = theta$mu, sigmas = theta$sigma, pie = theta$pie, z = z, 
              marginal_loglik = marginal.loglik)
}

# dgmm.loglik2 <- function (theta, z, marginal.loglik = FALSE) {  
#   TempFuncJoint <- function (k) {
#     theta$pie[k]*
#       dmvnormal(z, mu = theta$mu[[k]], sigma = theta$sigma[[k]])
#   }
#   loglik <- log(rowSums(rbind(sapply(1:theta$m, FUN = TempFuncJoint))))
#   if (!marginal.loglik)
#     loglik <- sum(loglik)
#   return(loglik)
# }