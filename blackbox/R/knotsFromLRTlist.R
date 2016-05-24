knotsFromLRTlist <- function(locLRTlist, lrthreshold) {
  ## adding key points to the training points: (redundant of) [the profile points for the LRT, the profile edges], and the global maximum
  keypts <- c() ## required void value in case there is no LRT point to be added
  if(length(locLRTlist)>0) {
    tmp <- t(as.data.frame(
      lapply(locLRTlist, function(l) {
        return(tofullKrigingspace(l$profpt$par, fixedlist=l$LRTfixedvals))
      })
    ))
    keyptscheck <- apply(tmp, 1, tofKpredict.nohull, fixedlist=NULL)
    keypts <- tmp[which(keyptscheck>lrthreshold), , drop=FALSE] ## in particular takes out LRT...$profpt$par with too low $value (necess for mirror->reflexion)
  }
  return(keypts)
}
