bounds.corr.GSC.PP <-
function (lamvec) 
{
  if (sum(lamvec <= 0) > 0) {
    stop("lambda should be positive \n")
  }
  norow = 1e+05
  maxmat = minmat = diag(NA, length(lamvec))
  errorCount = 0
  for (i in 2:length(lamvec)) {
    for (j in 1:(i-1)) { 
      Xpoisi = rpois(norow,lamvec[i])
      Xpoisj = rpois(norow,lamvec[j])
      max = cor(Xpoisi[order(Xpoisi)], Xpoisj[order(Xpoisj)])
      min = cor(Xpoisi[order(Xpoisi, decreasing = TRUE)], Xpoisj[order(Xpoisj)])
      minmat[i, j] = minmat[j, i] = min
      maxmat[i, j] = maxmat[j, i] = max
    }
  }
  return(list(min = round(minmat,3), max = round(maxmat,3)))
}
