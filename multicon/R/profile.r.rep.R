Profile.r.rep <-
function(x.set, y.set, nomiss=1.0, CI=.95, CItype="xci") {
  x.Norm <- colMeans(x.set, na.rm=T)
  y.Norm <- colMeans(y.set, na.rm=T)
  x.dist <- temp.resid(x.Norm, x.set, nomiss=nomiss)
  y.dist <- temp.resid(y.Norm, y.set, nomiss=nomiss)
  overall.rep <- alpha.cov(cov(ipsatize(x.set)*ipsatize(y.set), use="pair"))
  dist.rep <- alpha.cov(cov(ipsatize(x.dist)*ipsatize(y.dist), use="pair"))
  Reps <- c(overall.rep, dist.rep)
  if(CItype=="xci") {
    CIs <- sapply(Reps, alpha.xci, k=nrow(x.set), n=ncol(x.set), CI=CI)
  }
  if(CItype=="aci") {
    CIs <- sapply(Reps, alpha.aci, k=nrow(x.set), n=ncol(x.set), CI=CI)
  } 
  out <- t(rbind(Reps, CIs))
  colnames(out) <- c("Replicability", "Lower Limit", "Upper Limit")
  rownames(out) <- c("Overall", "Distinctive")
  return(out)
}
