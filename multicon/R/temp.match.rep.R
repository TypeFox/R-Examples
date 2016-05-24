temp.match.rep <-
function(template, y.set, CI=.95, CItype="xci") {
  y.Norm <- colMeans(y.set, na.rm=T)
  y.dist <- temp.resid(y.Norm, y.set, nomiss=1.0)
  temp.mat <- matrix(scale2(template), ncol=length(template), nrow=nrow(y.set), byrow=T)
  overall.rep <- alpha.cov(cov(ipsatize(y.set)*temp.mat, use="pair"))
  dist.rep <- alpha.cov(cov(ipsatize(y.dist)*temp.mat, use="pair"))
  Reps <- c(overall.rep, dist.rep)
  if(CItype=="xci") {
    CIs <- sapply(Reps, alpha.xci, k=nrow(y.set), n=ncol(y.set), CI=CI)
  }
  if(CItype=="aci") {
    CIs <- sapply(Reps, alpha.aci, k=nrow(y.set), n=ncol(y.set), CI=CI)
  } 
  out <- t(rbind(Reps, CIs))
  colnames(out) <- c("Replicability", "Lower Limit", "Upper Limit")
  rownames(out) <- c("Overall", "Distinctive")
  return(out)
}
