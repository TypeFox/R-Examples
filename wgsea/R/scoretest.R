scoretest <- function(x,y,n.perm=0,pheno.perm=NULL) {
  n <- length(y)
  if(!ncol(x) == n)
    stop("x should have dim p x n, where n is length(y).")

  if(n.perm>0)
    pheno.perm <- genperms(y,n.perm)
  if(!is.null(pheno.perm)) {
    if(nrow(pheno.perm) != n)
      stop("length of phenotype vector must equal total number of cases and controls.\n")
    n.perm <- ncol(pheno.perm)
  }
  
  permute <- n.perm > 0
  if(!permute) {
    return(pchisq(myscore(x,y),df=1,lower.tail=FALSE))
  }
  result <- matrix(NA,nrow(x),n.perm)
  for(i in 1:n.perm) {
    cat(".")
    result[,i] <- pchisq(myscore(x,pheno.perm[,i]),df=1)
  }
  return(result)
}
myscore <- function(x,y) { # x is a matrix of explanatory vars, y is a vector of response (0,1)
  x <- t(x - rowMeans(x))
  U <- colSums(y * x)
  n <- length(y)
  V <- sum( (y - mean(y))^2 ) * colSums(x^2) / n
  X <- U^2/V
  return(X)
}
