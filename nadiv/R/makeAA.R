makeAA <- function(pedigree)
{
  A <- makeA(pedigree)
  AA <- A*A
  logDet <- determinant(AA, logarithm = TRUE)$modulus[1]
  AAinv <- solve(AA)
  listAAinv <- sm2list(AAinv, rownames=pedigree[,1], colnames=c("row", "column", "AAinverse"))
 return(list(AA = AA, logDet = logDet, AAinv = AAinv, listAAinv = listAAinv))    
}

