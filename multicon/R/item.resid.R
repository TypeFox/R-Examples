item.resid <-
function(x.set, y.set, nomiss=.8) {
  Ns <- apply(x.set*y.set, 2, function(x) length(x) - sum(is.na(x)))
  b1 <- (colSums(scale2(x.set)*scale2(y.set)) / Ns) * (apply(y.set, 2, sd, na.rm=T) / apply(x.set, 2, sd, na.rm=T))
  b0 <- colMeans(y.set, na.rm=TRUE) - b1*colMeans(x.set, na.rm=TRUE)
  b1 <- ifelse(Ns < nomiss*nrow(y.set), NA, b1)
  b0 <- ifelse(Ns < nomiss*nrow(x.set), NA, b0)
  b1.mat <- matrix(b1, nrow=nrow(y.set), ncol=ncol(y.set), byrow=TRUE)
  b0.mat <- matrix(b0, nrow=nrow(y.set), ncol=ncol(y.set), byrow=TRUE)
  pred.mat <- b0.mat + b1.mat*x.set
  res <- data.frame(y.set - pred.mat)
  return(res)
}
