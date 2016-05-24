detvarorglm <- function(object){
  eta <- object$X %*% object$coefficients
  y <- object$y
  z <- as.vector(eta  + (y - object$family$linkinv(eta))/object$family$mu.eta(eta))
  w <- sqrt(object$weights)  
  xm <- object$X*w
  A <- object$constr
  lmf <- lm(z * w ~ xm-1)  
  
  sx <- summary(lmf)$cov.unscaled
  sA <- qr.solve(A %*% sx %*% t(A))
  M <- diag(ncol(xm)) - sx %*% t(A) %*% sA %*% A
  Ms <- (diag(ncol(xm)) - t(A) %*% sA %*% A %*% sx)
  vc <- M %*% sx %*% Ms
  
  svdd <- svd(vc)$d
  eok <- svdd > 2 * .Machine$double.eps  
  
  evv <- svd(vc[A!=0,A!=0])$d
  evok <- evv > 2 * .Machine$double.eps  
  
  dvv <- prod(evv[evok])
  dvca <- prod(svdd[eok])/dvv
  return(dvca)
}

