weight.boots <- function(object, ncomp = object$ncomp, conf = .95) {
  if ((object$val.method == "none" | object$val.method == "loo")) {
    stop("No bootstrapping was done for this model")
  }
  conf <- conf
  weights.boot.a <- do.call("rbind", object$validation$weights)
  weights.boot <- as.matrix(weights.boot.a)
  Upper <- 1 - (((1 - conf)/2))
  Lower <- 1 - Upper
  weights.boot.cis <- llply(ncomp, function(y) {
    do.call("rbind", as.list(
      by(weights.boot[, y], list(row.names(weights.boot)), function(x){
        c(ncomp = y, boot.mean = mean(x, na.rm = T), Skewness = skewness(x, na.rm = T), 
          quantile(x, c(Lower, Upper), na.rm = T), 'Bootstrap Error' = sd(x, na.rm = T))
      }
      )))
  })
  names(weights.boot.cis) <- ncomp
  OC <- data.frame(object$weights[, ncomp])
  names(OC) <- ncomp
  A <- llply(1:length(weights.boot.cis), function(x) {
    x. <- names(weights.boot.cis)[x]
    weights.boot.cis2 <- as.data.frame(weights.boot.cis[[x]])
    weights.boot.cis2$variable <- row.names(weights.boot.cis[[x]])
    row.names(weights.boot.cis2) <- NULL
    Out <- weights.boot.cis2[as.factor(row.names(object$weights)), ]
    Out$Actual <- OC[, x.]
    Out$Actual <- object$weights[, x]
    Out$Bias <- Out$boot.mean - Out$Actual
    Out$'t value' <- Out$Actual / Out$'Bootstrap Error'
    Out$'bias t value' <- Out$Bias / Out$'Bootstrap Error'
    Out[, c(1, 7, 8, 4:5, 2, 3, 9, 6, 10:11)]
  })
  A
}