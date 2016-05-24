Profile.r <-
function(x.set, y.set, nomiss=1.0, distinct=FALSE, alt="greater") {
  if(distinct==F) {
    valids <- apply(cbind(x.set, y.set), 1, function(x) valid.pairs(x[1:ncol(x.set)], x[-1:-ncol(y.set)])$Pct)
    out <- diag(cor(t(x.set),t(y.set), use="pair"))
    out <- ifelse(valids >= nomiss, out, NA)
  }
  if(distinct==T) {
    valids <- apply(cbind(x.set, y.set), 1, function(x) valid.pairs(x[1:ncol(x.set)], x[-1:-ncol(y.set)])$Pct)
    mat <- cor(t(x.set),t(y.set), use="pair")
    overall.r <- diag(mat)
    overall.r <- ifelse(valids >= nomiss, overall.r, NA)
    baseline <- mean(fisherz(mat[upper.tri(mat)]), na.rm=TRUE)
    x.Norm <- colMeans(x.set, na.rm=T)
    y.Norm <- colMeans(y.set, na.rm=T)
    norm.r <- cor(x.Norm, y.Norm)
    x.dist <- temp.resid(x.Norm, x.set, nomiss=nomiss)
    y.dist <- temp.resid(y.Norm, y.set, nomiss=nomiss)
    dist.r <- diag(cor(t(x.dist),t(y.dist), use="pair"))
    dist.r <- ifelse(valids >= nomiss, dist.r, NA)
    overall.test <- t.test(fisherz(overall.r), mu=baseline, alternative=alt)
    dist.test <- t.test(fisherz(dist.r), alternative=alt)
    overall.T <- rbind(length(overall.r) - sum(is.na(overall.r)), fisherz2r(mean(fisherz(overall.r), na.rm=T)), fisherz2r(baseline), overall.test$statistic, overall.test$p.value)
    dist.T <- rbind(length(dist.r) - sum(is.na(dist.r)), fisherz2r(mean(fisherz(dist.r), na.rm=T)), 0, dist.test$statistic, dist.test$p.value)
    tests <- data.frame(overall.T, dist.T, row.names=c("N", "Mean", "baseline", "t", "p-value"))
    colnames(tests) <- c("Overall", "Distinctive")
    agrees <- data.frame("Overall"=overall.r, "Distinctive"=dist.r, row.names=rownames(x.set))
    out <- list("xNorm"=x.Norm, "yNorm"=y.Norm, "Norm.r"=norm.r, "Agreement"=agrees, "Tests"=tests)
  }
  return(out)
}
