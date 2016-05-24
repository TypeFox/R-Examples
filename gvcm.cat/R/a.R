a <- function (p = NULL, m = 1, ...)
{ k <- length(p)
  m <- min(k-1,m)
  if (k>1)
  {
    if (k > 2)
    {
      mat <- matrix(nrow=k, ncol=0)
      for (j in 1:m){
        mat2 <- diag(1, ncol=k, nrow=k)
        for (i in 1:(k-j)){
          mat2[p[i],p[i+j]] <- -1
        }
        mat2 <- mat2[,-p[1:j]]
        mat <- cbind (mat, mat2)
      }
    } else { mat <- matrix(c(-1,1), ncol=1) }  # k == 2
  } else {
    mat <- matrix(ncol=0,nrow=1)
  }
  colnames(mat) <- c()
  return (mat)
}
