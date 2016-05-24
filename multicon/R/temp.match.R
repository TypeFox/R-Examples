temp.match <-
function(template, y.set, nomiss=1.0, distinct=FALSE) {
  valids <- apply(y.set, 1, function(x) (length(x) - sum(is.na(x))) / length(x))
  if(distinct==F) {
    out <- as.vector(cor(template, t(y.set), use="pair"))
    out <- ifelse(valids >= nomiss, out, NA)
  }
  if(distinct==T) {
    overall.r <- as.vector(cor(template, t(y.set), use="pair"))
    overall.r <- ifelse(valids >= nomiss, overall.r, NA)
    y.Norm <- colMeans(y.set, na.rm=T)
    y.dist <- temp.resid(y.Norm, y.set, nomiss=nomiss)
    dist.r <- cor(template, t(y.dist), use="pair")
    dist.r <- ifelse(valids >= nomiss, dist.r, NA)
    matches <- data.frame("Overall"=overall.r, "Distinctive"=dist.r, row.names=rownames(y.set))
    out <- list("yNorm"=y.Norm, "Matches"=matches)
  }
  return(out)
}
