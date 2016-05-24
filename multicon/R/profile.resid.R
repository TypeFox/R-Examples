Profile.resid <-
function(x.set, y.set, nomiss=.8) {
  coefs <- Profile.reg(x.set=x.set, y.set=y.set, center="none", nomiss=nomiss)
  res <- data.frame(y.set - (coefs[,1] + coefs[,2]*x.set))
  names(res) <- paste(colnames(y.set), sep="", ".res")
  return(res)
}
