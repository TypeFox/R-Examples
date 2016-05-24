## Tools for manipulation of text count matrices ##

## converting count to frequency matrix
normalize <- function(x, byrow=TRUE){
    if(byrow){ s <- row_sums(x)
               s[s==0] <- 1
               return( x/s ) }
    else{
      s <- col_sums(x)
      s[s==0] <- 1
      return(t(t(x)/s)) }
}

## converting a count/freq matrix to tfidf
stm_tfidf <- function(x){
  idf <- log( nrow(x) ) - log(col_sums(x>0) + 1) 
  t( t(x) * idf )
}
    
## Dirichlet RNG
rdir <- function(n, alpha)
{
    x <- matrix(rgamma(length(alpha)*n,alpha),nrow=n,byrow=TRUE)
    return(t(x/rowSums(x))) }
