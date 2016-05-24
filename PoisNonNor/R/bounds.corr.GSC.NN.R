bounds.corr.GSC.NN <-
function (pmat) 
{
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
  
  maxmat = minmat = diag(NA, dim(pmat)[1])
  for (i in 2:dim(pmat)[1]) {
    for (j in 1:(i-1)) { 
      Yi = fleishman.uni(matrix(pmat[i,],nrow=1))
      Yj = fleishman.uni(matrix(pmat[j,],nrow=1))
      max = cor(Yi[order(Yi)], Yj[order(Yj)])
      min = cor(Yi[order(Yi, decreasing = TRUE)], Yj[order(Yj)])
      minmat[i, j] = minmat[j, i] = min
      maxmat[i, j] = maxmat[j, i] = max
    }
  }
  return(list(min = round(minmat,3), max = round(maxmat,3)))
}
