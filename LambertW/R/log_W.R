#' @rdname W
#' @examples
#' curve(log(x), 0.1, 5, lty = 2, col = 1, ylab = "")
#' curve(W(x), 0, 5, add = TRUE, col = "red")
#' curve(log_W(x), 0.1, 5, add = TRUE, col = "blue")
#' grid()
#' legend("bottomright", c("log(x)", "W(x)", "log(W(x))"),
#'        col = c("black", "red", "blue"), lty = c(2, 1, 1))
#' 
#' @export
log_W <- function(z, branch = 0, W.z = W(z, branch = branch)) {
  # log(W(z)) = log(z) - W(z)
  log.W.z <- log(z) - W.z
  if (any(is.infinite(z))) {
    log.W.z[is.infinite(z)] <- Inf
  }
  return(log.W.z)
} 

