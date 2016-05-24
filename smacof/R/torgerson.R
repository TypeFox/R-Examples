torgerson <- function(delta, p = 2)
{
  
#diss ... dissimilarity matrix
#p ... number of dimensions
#------------------- begin subroutines -------------------
#Torgerson's double centering
  doubleCenter <- function(x) {
    n <- dim(x)[1]
    m <- dim(x)[2]
    s <- sum(x)/(n*m)
    xr <- rowSums(x)/m
    xc <- colSums(x)/n
    return((x-outer(xr,xc,"+"))+s)
  }
#-------------------- end subroutines --------------------
  diss <- as.matrix(delta)
  
  z <- eigen(-doubleCenter(as.matrix(diss)^2)/2,symmetric=TRUE)
  v <- pmax(z$values,0)
  if (p == 1) normdiag <- cbind(sqrt(v[1])) else normdiag <- diag(sqrt(v[1:p]))
  conf <- z$vectors[,1:p]%*%normdiag
  rownames(conf) <- rownames(diss)
  return(conf)
}

