coefficients.boots <- function(object, ncomp = object$ncomp, conf = .95) {
  if ((object$val.method == "none" | object$val.method == "loo")) {
    stop("No bootstrapping was done for this model")
  }
  conf <- conf
  coefficients.boot.a <- do.call("rbind", object$validation$coefficients)
  coefficients.boot <- as.matrix(coefficients.boot.a)
  Upper <- 1 - (((1 - conf)/2))
  Lower <- 1 - Upper
  coefficients.boot.cis <- llply(ncomp, function(y) {
    do.call("rbind", as.list(
      by(coefficients.boot[, y], list(row.names(coefficients.boot)), function(x){
        c(ncomp = y, boot.mean = mean(x, na.rm = T), Skewness = skewness(x, na.rm = T), 
          quantile(x, c(Lower, Upper), na.rm = T), 'Bootstrap Error' = sd(x, na.rm = T))
      }
      )))
  })
  names(coefficients.boot.cis) <- ncomp
  OC <- data.frame(object$coefficients[, ncomp])
  names(OC) <- ncomp
A <- llply(1:length(coefficients.boot.cis), function(x) {
    x. <- names(coefficients.boot.cis)[x]
    coefficients.boot.cis2 <- as.data.frame(coefficients.boot.cis[[x.]])
    coefficients.boot.cis2$variable <- row.names(coefficients.boot.cis[[x.]])
    row.names(coefficients.boot.cis2) <- NULL
    Out <- coefficients.boot.cis2[as.factor(row.names(OC)), ]
    Out$Actual <- OC[, x.]
    Out$Bias <- Out$boot.mean - Out$Actual
    Out$'t value' <- Out$Actual / Out$'Bootstrap Error'
    Out$'bias t value' <- Out$Bias / Out$'Bootstrap Error'
    Out[, c(1, 7, 8, 4:5, 2, 3, 9, 6, 10:11)]
  })
  A
}
