R2s <- function(object) {
  if(is.null(object$validation)) {
    stop("For this to work you need to cross-validate; set 'validation' = 'oob' or 'loo'")
  }
  R2X <- sapply(1:object$ncomp, function(x) {
    1 - sum(diag(crossprod(as.matrix(object$Xdata) - 
                             (object$scores[ , 1:x] %*% t(object$loadings[ , 1:x]))))) / 
      sum(diag(crossprod(as.matrix(object$Xdata))))
  })
  R2Y <- 1 - sapply(1:object$ncomp, function(x) crossprod(object$Yactual - object$iPreds[, x])) / 
    crossprod(object$Yactual - mean(object$Yactual))
  R2.values <- list(Comp = 1:length(as.vector(object$validation$cvR2)), 
                          CVR2 = as.vector(object$validation$cvR2), 
                          R2X = R2X, R2Y = R2Y)
  class(R2.values) <- "R2s"
  R2.values
}

print.R2s <- function(x, ncomp = x$VIP$ncomp, ...) {
  R2s.output <- data.frame(Comp = x$Comp, 
             CVR2 = x$CVR2, 
             R2X = x$R2X, R2Y = x$R2Y)
  print(R2s.output)
}

