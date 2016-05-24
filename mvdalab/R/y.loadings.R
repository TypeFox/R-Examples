y.loadings <- y.loadings.boots <- function (object, conf = .95) {
  if (object$val.method == "none" | object$val.method == "loo") {
    B <- object$y.loadings
    return(B)
  } else {
    conf <- conf
    y.loadings.boots <- do.call("rbind", object$validation$y.loadings)
    Upper <- 1 - (((1 - conf)/2))
    Lower <- 1 - Upper
    y.loadings.boot.cis <- t(apply(y.loadings.boots, 1, function(x) quantile(x, c(Lower, Upper), na.rm = T)))
    row.names(y.loadings.boot.cis) <- paste("ncomp", 1:object$ncomp, sep = " ")
    Out <- data.frame(Actual = object$y.loadings[, 1], y.loadings.boot.cis)
    names(Out)[2:3] <- paste(c(Lower * 100, Upper * 100), "%", sep = "")
    Out$boot.mean <- apply(y.loadings.boots, 1, function(x) mean(x, na.rm = T))
    Out$Skewness <- apply(y.loadings.boots, 1, function(x) skewness(x, na.rm = T))
    Out$Bias <- Out$boot.mean - Out$Actual
    Out$'Bootstrap Error' <- apply(y.loadings.boots, 1, function(x) sd(x, na.rm = T))
    Out$'t value' <- Out$Actual / Out$'Bootstrap Error'
    Out$'bias t value' <- Out$Bias / Out$'Bootstrap Error'
    return(Out)
  }
}