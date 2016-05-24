#' @rdname W
#' @export

deriv_W <- function(z, branch = 0, W.z = W(z, branch = branch)) {
  
  if (branch == 0) {
    W.deriv <- exp(log_deriv_W(z, branch = branch, W.z = W.z))
    # W.deriv <- exp(log_deriv- W.z - log(W.z + 1))
  } else if (branch == -1) {
    # for negative branch, W(z, -1) << -1 is possible, hence log(W.z + 1) 
    # would be zero.  Hence we can't use the log transform to simplify 
    # expression
    W.deriv <- rep(NA, length(z))
    W.deriv[z == 0] <- 1
    ind.neg <- (z < 0)
    W.deriv[ind.neg] <- (1 - 1/(W.z[ind.neg] + 1)) / z[ind.neg]
    dim(W.deriv) <- dim(z)
  }
  return(W.deriv)
} 

#' @rdname W
#' @export
#' 

log_deriv_W <- function(z, branch = 0, W.z = W(z, branch = branch)) {
 
  log.deriv.W <- rep(NA, length(z))
  ind.neg <- z < 0
  # for z >= 0 there is a closed form
  log.deriv.W[!ind.neg] <- - W.z[!ind.neg] - log(W.z[!ind.neg] + 1)
  # for z < 0, the closed form also exists but since both are negative the
  # log() would return NA.  Use a trick by multiplying both expressions by -1.
  # W'(z) = (1 - 1 / (W(z) + 1)) * z
  # For logarithms this can be simplified by multiplying both by -1
  if (branch == 0) {
    log.deriv.W[ind.neg] <- log(-(1 - 1/(W.z[ind.neg] + 1))) - log(-z[ind.neg])
  } else if (branch == -1) {
    log.deriv.W[ind.neg] <- log(deriv_W(z = z[ind.neg], 
                                        branch = branch, W.z = W.z[ind.neg]))
  } else {
    stop("Branch must be either 0 or -1.")
  }
  return(log.deriv.W)
}


#' @rdname W
#' @export
#' 

deriv_log_W <- function(z, branch = 0, W.z = W(z, branch = branch)) {
    
  deriv.log.W <-  (z * (1 + W.z))^-1
  deriv.log.W[z < 0] <- NaN
  return(deriv.log.W)
}
