#' @rdname W_delta
#' @export
deriv_W_delta <- function(z, delta = 0) {
  
  stopifnot(is.numeric(z), 
            is.numeric(delta), 
            length(delta) == 1)
  
  out <- rep(1, length(z))
  if (delta != 0) {
    ind.zero <- (z == 0)
    value <- 1
    value[!ind.zero] <- W(delta * z[!ind.zero]^2)
    out[!ind.zero] <- 
      delta^(-0.5) * value[!ind.zero]^(-0.5) * 
        value[!ind.zero]/(z[!ind.zero] * (1 + value[!ind.zero]))
    out[!ind.zero] <- sign(z[!ind.zero]) * out[!ind.zero]
  }
  return(out)
} 

#' @rdname W_delta
#' @export
deriv_W_delta_alpha <- function(z, delta = 1, alpha = 1) {
  
  stopifnot(is.numeric(z),
            is.numeric(delta),
            is.numeric(alpha),
            length(delta) == 1,
            length(alpha) == 1)
  
  out <- rep(1, length(z))
  if (delta != 0) {
    da <- delta * alpha
    W.da <- W(da * (z^2)^alpha)
    
    first.term <- (W.da/da)^(1/(2 * alpha) - 1)
    second.term <- W.da/(da * z * (1 + W.da))
    
    out <- sign(z) * first.term * second.term
    out[is.na(out)] <- 1
  } 
  return(out)
}