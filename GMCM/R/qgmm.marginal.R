#' @rdname dgmcm.loglik
qgmm.marginal <- function (u, theta, res = 1000, spread = 5) {
  d <- theta$d
  m <- theta$m
  n.samples <- round(res * theta$pie)
  n.samples[n.samples == 0] <- 2
  
  # Create grid of evalutation
  s <- NULL
  for (i in 1:d) {
    for (j in 1:m) {
      m.ij <- theta$mu[[j]][i]
      sd.ij <- sqrt(theta$sigma[[j]][i,i])
      s <- c(s, seq(m.ij-spread*sd.ij, m.ij+spread*sd.ij, l = n.samples[j]))
    }
  }
  dim(s) <- c(sum(n.samples), d)
  
  # Evaluate on cdf on the grid
  eval <- 
    pgmm_marginal(z = s, mus = theta$mu, sigmas = theta$sigma, pie = theta$pie)
  
  # Invert function
  z.out <- NULL
  for (j in 1:d) {
    z.out <- c(z.out, approxfun(eval[, j], s[, j], rule = 2)(u[, j]))
  }
  z.out.is.na <- is.na(z.out)
  if (any(z.out.is.na)) {
    z.out[z.out.is.na & u >= 1] <- Inf
    z.out[z.out.is.na & u <= 0] <- -Inf
  }
  dim(z.out) <- c(nrow(u),d)
  return(z.out)
}

# qgmm.marginal2 <- function (u, theta, res = 1000, spread = 5) 
# {
#   
#   d <- theta$d
#   m <- theta$m
#   n.samples <- round(res * theta$pie)
#   n.samples[n.samples == 0] <- 2
#   
#   # Create grid of evalutation
#   s <- NULL
#   for (i in 1:d) {
#     for (j in 1:m) {
#       m.ij <- theta$mu[[j]][i]
#       sd.ij <- sqrt(theta$sigma[[j]][i,i])
#       s <- c(s, seq(m.ij-spread*sd.ij, m.ij+spread*sd.ij, l = n.samples[j]))
#     }
#   }
#   dim(s) <- c(sum(n.samples), d)
# 
#   # Evaluate on cdf on the grid
#   eval <- pgmm.marginal(z = s, theta = theta)
#   
#   # Invert function
#   z.out <- NULL
#   for (j in 1:d) {
#     z.out <- c(z.out, approxfun(eval[, j], s[, j], rule = 2)(u[, j]))
#   }
#   if (any(is.na(z.out))) {
#     z.out[is.na(z.out) & u >= 1] <- Inf
#     z.out[is.na(z.out) & u <= 0] <- -Inf
#   }
#   dim(z.out) <- c(nrow(u),d)
#   return(z.out)
# }
