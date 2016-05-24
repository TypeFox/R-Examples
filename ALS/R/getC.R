"getC" <- function(EAllAllTimes, PsiAll_T, C, baseline, uni, nonnegC, x) { 
  for(i in 1:ncol(PsiAll_T)) 
    C[i,] <- if(nonnegC) coef(nnls::nnls(A = EAllAllTimes[[i]], b = PsiAll_T[,i]))
    else qr.coef(qr( EAllAllTimes[[i]]), PsiAll_T[,i])
  if(uni) {
    ncolel <- ncol(C)
    if(baseline)
      ncolel <- ncolel - 1 
    for(i in 1:ncolel) {
      C[,i] <- Iso::ufit(y=C[,i],x=x)$y
    }
  }
  C
}

