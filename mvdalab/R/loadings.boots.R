loadings.boots <- function(object, ncomp = object$ncomp, conf = .95) {
  if ((object$val.method == "none" | object$val.method == "loo")) {
    stop("No bootstrapping was done for this model")
  }
  conf <- conf
  loadings.boot.a <- do.call("rbind", object$validation$loadings)
  loadings.boot <- as.matrix(loadings.boot.a)
  Upper <- 1 - (((1 - conf)/2))
  Lower <- 1 - Upper
  loadings.boot.cis <- llply(ncomp, function(y) {
    do.call("rbind", as.list(
      by(loadings.boot[, y], list(row.names(loadings.boot)), function(x){
        c(ncomp = y, boot.mean = mean(x, na.rm = T), Skewness = skewness(x, na.rm = T), 
          quantile(x, c(Lower, Upper), na.rm = T), 'Bootstrap Error' = sd(x, na.rm = T))
      }
      )))
  })
  names(loadings.boot.cis) <- ncomp
  OC <- data.frame(object$loadings[, ncomp])
  names(OC) <- ncomp
A <- llply(1:length(loadings.boot.cis), function(x) {
    x. <- names(loadings.boot.cis)[x]
    loadings.boot.cis2 <- as.data.frame(loadings.boot.cis[[x.]])
    loadings.boot.cis2$variable <- row.names(loadings.boot.cis[[x.]])
    row.names(loadings.boot.cis2) <- NULL
    Out <- loadings.boot.cis2[as.factor(row.names(object$loadings)), ]
    Out$Actual <- OC[, x.]
    Out$Actual <- object$loadings[, x]
    Out$Bias <- Out$boot.mean - Out$Actual
    Out$'t value' <- Out$Actual / Out$'Bootstrap Error'
    Out$'bias t value' <- Out$Bias / Out$'Bootstrap Error'
    Out[, c(1, 7, 8, 4:5, 2, 3, 9, 6, 10:11)]
  })
  A
}


