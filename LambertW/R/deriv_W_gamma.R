#' @rdname W_gamma
#' @export
deriv_W_gamma <- function(z, gamma = 0, branch = 0) {
  stopifnot(is.numeric(z), 
            is.numeric(gamma), 
            length(gamma) == 1)
  
  out <- rep(1, length(z))
  if (gamma != 0) {
    out <- deriv_W(gamma * z, branch = branch)
  }
  dim(out) <- dim(z)
  return(out)
} 
