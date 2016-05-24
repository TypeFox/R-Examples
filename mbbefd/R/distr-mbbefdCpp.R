#' @useDynLib mbbefd
#' @importFrom Rcpp sourceCpp


#random generation
rmbbefdCpp <- function(n, a, b)
{
  .rmbbefdC(n,  a,  b) 
}

rMBBEFDCpp <- function(n, g, b)
{
  .rMBBEFDC(n, g, b) 
}
  

### r function MBBEFD(a,b) for users
rmbbefd <- rmbbefdCpp

### r function MBBEFD(g,b)
rMBBEFD <- rMBBEFDCpp

