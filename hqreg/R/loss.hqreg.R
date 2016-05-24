loss.hqreg <- function(y, yhat, args) {
  r <- y-yhat
  type.measure <- args$type.measure
  if (type.measure == "deviance") {
    method <- args$method
    if (method == "huber") {
      gamma <- args$gamma
      val <- hloss(r, gamma)
    } else if (method == "quantile") {
      tau <- args$tau
      val <- qloss(r, tau)
    } else {
      val <- r^2
    }    
  } else if (type.measure == "mse") {
    val <- r^2
  } else {
    val <- abs(r)
  }
  val
}

hloss <- function(r, gamma) {
  rr <- abs(r)
  ifelse(rr <= gamma, rr^2/(2*gamma), rr-gamma/2)
}

qloss <- function(r, tau) ifelse(r <= 0, (tau-1)*r, tau*r)
