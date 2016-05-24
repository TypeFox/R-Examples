"getS" <- function(CList, PsiAll, S, W, baseline, uni, nonnegS, normS, x2) {
  C <- do.call("rbind",CList)
  C[which(is.nan(C))] <- 1
  for(i in 1:ncol(PsiAll)) {
    if(nonnegS)
      s <- try(nnls::nnls(A = C * W[,i], b = PsiAll[,i]))
    else
      s <- try(qr.coef( qr(C * W[,i]), PsiAll[,i]))
    if(class(s) == "try-error")
      S[i,] <- rep(1, ncol(C))
    else S[i,] <- if(nonnegS) coef(s) else s
  }
  if(uni) {
    ncolel <- ncol(C)
    if(baseline)
      ncolel <- ncolel - 1 
    for(i in 1:ncolel) 
      S[i,] <- Iso::ufit(y=S[i,],x=x2)$y
  }
  if(normS>0) {
    if(normS==1) 
      S <- normdat(S)
    else {
      for(i in 1:ncol(S)) {
        nm <- sqrt(sum((S[,i])^2))
        S[,i] <- S[,i]/nm
      }
    }         
  }
  S
}

