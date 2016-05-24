Profile.norm <-
function(set, nomiss=.8, center="grand", std=FALSE) {
  colSum <- colSums(set, na.rm=T)
  colLen <- apply(set, 2, length)
  colMiss <- apply(set, 2, function(x) sum(is.na(x)))
  colVal <- colLen - colMiss
  sum.mat <- matrix(colSum, nrow=nrow(set), ncol=ncol(set), byrow=T)
  jackMeans <- (sum.mat - set) / (colVal - 1)
  cors <- Profile.r(set, jackMeans, nomiss=nomiss)
  regs <- Profile.reg(x.set=jackMeans, y.set=set, nomiss=nomiss, center=center, std=std)
  resid.df <- Profile.resid(x.set=jackMeans, y.set=set, nomiss=nomiss)
  out <- list("Means"=colMeans(set, na.rm=T), "JackMeans"=jackMeans, "Cors"=cors, "Regs"=regs, "Residuals"=resid.df)
  return(out)
}
