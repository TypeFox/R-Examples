"getCList" <- function(S, PsiList, CList, WList, resid, x, baseline,
                       fixed, uni, nonnegC, closureC) {
  for(j in 1:length(PsiList)) {
    S[which(is.nan(S))] <- 1
    if(length(fixed[[j]])>0) 
      S <- S[, -fixed[[j]]]
    for(i in 1:nrow(PsiList[[j]])) {
      if(nonnegC)
        cc <- try(nnls::nnls(A = S * WList[[j]][i,], b = PsiList[[j]][i,]))
      else
        cc <- try(qr.coef(qr(S * WList[[j]][i,]), PsiList[[j]][i,]))
      if(class(cc) == "try-error")
        sol <- rep(1, ncol(S))
      else
        sol <- if(nonnegC) coef(cc) else cc 
      cc1 <- rep(NA, ncol(CList[[j]]))  
      if(length(fixed[[j]])>0) 
        cc1[fixed[[j]]] <- 0
      cc1[is.na(cc1)] <- sol 
      CList[[j]][i,] <- cc1 
    }
  }
  if(uni) {
    for(j in 1:length(PsiList)) {
      ncolel <- ncol(CList[[j]])
      if(baseline)
        ncolel <- ncolel - 1 
      for(i in 1:ncolel) {
        CList[[j]][,i] <- Iso::ufit(y=CList[[j]][,i],x=x)$y
      }
    }
  }
  if(length(closureC) > 1) {
    for(j in 1:length(PsiList)) 
      for(i in 1:nrow(PsiList[[j]]))
        CList[[j]][i,] <- sum((CList[[j]][i,]*closureC[[j]][i])/
                              max(sum(CList[[j]][i,])))
  }
  CList
}

