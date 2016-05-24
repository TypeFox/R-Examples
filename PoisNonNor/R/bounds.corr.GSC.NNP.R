bounds.corr.GSC.NNP <-
function (lamvec, pmat) {
  if (sum(lamvec <= 0) > 0) {
    stop("lambda should be positive \n")
  }

  if (dim(pmat)[2] != 4){
    stop("column of pmat must be 4\n")
  }
  
  fleishman.uni = function(p,norow=1e+05)
  {
    x = rnorm(norow)
    X = as.matrix(cbind(rep(1,norow),x,x^2,x^3))
    Y = X%*%t(p)
    return(Y)
  }
  
  norow = 1e+05
  maxmat = minmat = matrix(NA, nrow=length(lamvec),ncol=dim(pmat)[1])
  errorCount = 0
  for (i in 1:length(lamvec)) {
    for (j in 1:dim(pmat)[1]) { 
      Xpoisi = rpois(norow,lamvec[i])
      Yj = fleishman.uni(matrix(pmat[j,],nrow=1))
      max = cor(Xpoisi[order(Xpoisi)], Yj[order(Yj)])
      min = cor(Xpoisi[order(Xpoisi, decreasing = TRUE)], Yj[order(Yj)])
      minmat[i, j] = min
      maxmat[i, j] = max
    }
  }
  return(list(min = round(minmat,3), max = round(maxmat,3)))
}
