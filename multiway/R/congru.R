congru <- function(x, y = NULL){
  ###### Tucker's Congruence Coefficient
  ###### Nathaniel E. Helwig (helwig@umn.edu)
  ###### Last modified: May 10, 2015
  
  x <- as.matrix(x)
  if(is.null(y)){
    if(ncol(x)==1L){
      stop("supply both 'x' and 'y' or a matrix-like 'x'")
    } else {
      dx <- diag(1/sqrt(colSums(x^2)))
      return( dx %*% crossprod(x) %*% dx )
    }
  } else {
    y <- as.matrix(y)
    if(nrow(x) != nrow(y)){stop("inputs 'x' and 'y' must have same dimensions")}
    ncx <- ncol(x)
    ncy <- ncol(y)
    if(ncx==1L & ncy==1L){
      return( sum(x*y) / sqrt( sum(x^2)*sum(y^2) ) )
    } else {
      if(ncx==1L) { dx <- matrix(1/sqrt(sum(x^2))) } else { dx <- diag(1/sqrt(colSums(x^2))) }
      if(ncy==1L) { dy <- matrix(1/sqrt(sum(y^2))) } else { dy <- diag(1/sqrt(colSums(y^2))) }
      return( dx %*% crossprod(x,y) %*% dy )
    }
  }
  
}