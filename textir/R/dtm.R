## Tools for manipulation of text count matrices ##
stm2dg <- function(x, ...){
  sparseMatrix(i=x$i, j=x$j, x=x$v,
        dims=dim(x),dimnames=dimnames(x), ...) }

## correlation 
corr <- function(x, y){
  if(inherits(x, "simple_triplet_matrix")) x <- stm2dg(x)
  if(!inherits(x, "dgCMatrix")){ return(cor(x,y) ) }

  n <- nrow(x)
  v <- t(scale(y))
  
  r <- tcrossprod(t(x)/sdev(x), v)/(n-1)
  dimnames(r) <- list(dimnames(x)[[2]], dimnames(y)[[2]])
  return( r ) 
}
  
## column standard deviation 
sdev <- function(x){
  if(is.null(dim(x))) return(sd(x))
  if(inherits(x, "simple_triplet_matrix")) x <- stm2dg(x)
  if(!inherits(x, "dgCMatrix")){ return(apply(as.matrix(x),2,sd)) }
  n <- nrow(x)
  s <- sqrt(colSums(x^2)/(n-1) - colSums(x)^2/(n^2 - n)) 
  names(s) <- colnames(s)
  return(s)  }

## tfidf
tfidf <- function(x, normalize=TRUE){
   if(inherits(x, "simple_triplet_matrix")) x <- stm2dg(x)      
   idf <- log(nrow(x)) - log(colSums(x>0) + 1)
   if(normalize) x <- x/rowSums(x)
   t( t(x)*idf )
}