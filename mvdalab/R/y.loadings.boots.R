y.loadings.boots <- function(object, ncomp = object$ncomp, conf = .95) {
  if ((object$val.method == "none" | object$val.method == "loo")) {
    stop("No bootstrapping was done for this model")
  }
  conf <- conf
  Q.boots <- as.matrix(object$validation$y.loadings[, 1:ncomp])
  Q.boot.cis <- lapply(1:ncol(Q.boots), function(x) {
    Upper <- 1 - (((1 - conf)/2))
    Lower <- 1 - Upper
    Percentiles <- quantile(Q.boots[, x], c(Lower, Upper), na.rm = T)
    boot.results <- Percentiles
    boot.results
  })
  do.call("rbind", Q.boot.cis)
}



