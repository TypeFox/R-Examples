horn <-
function(set, sims=100, nomiss=1.0, graph=TRUE) {
  Miss <- apply(set, 1, function(x) sum(is.na(x)))
  PercValid <- (ncol(set) - Miss) / ncol(set)
  comp <- set[PercValid >= nomiss,]
  cases.del <- nrow(set) - nrow(comp)
  value.out <- matrix(NA, nrow=sims, ncol=ncol(comp))
  for(i in 1:sims) {
    mat <- matrix(sample(unlist(comp), size=nrow(comp)*ncol(comp), replace=T), nrow=nrow(comp), ncol=ncol(comp))
    value.out[i,] <- eigen(cor(mat, use="pair"))$values
  }
  n.comps <- which(!(eigen(cor(comp, use="pair"))$values > colMeans(value.out)))[1] - 1
  if(graph==T) {
    op <- par(las=1, font.main=1)
    plot(eigen(cor(comp, use="pair"))$values, type="b", main="Scree Plot for Parallel Analysis", xlab="Number of Components", ylab="Eigenvalue")
    lines(colMeans(value.out), type="l", col="red")
    legend("topright", c("Simulated Values"), col="red", lty=1, bty="n")
  }
  cat("Parallel analysis suggests", n.comps, "components \n")
  cat(cases.del, "cases deleted due to missingness. \n")
}
