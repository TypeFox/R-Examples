#' Randomly Generate Low-Rank Correlation Matrix
#' 
#' Generate a correlation matrix as R = LL' where the rows of L are of length
#' 1, L is of rank \code{r} and the matrix L is sparse (depending on
#' \code{sparse.prop}. The loadings in L are sampled from a standard normal
#' distribution, after which \code{sparse.prop} is used to set a randomly
#' chosen number of loadings in each row equal to zero. To ensure that a
#' correlation matrix results, the rows are normalized.
#' 
#' @param m integer; the number of variables.
#' @param r integer; the required rank.
#' @param sparse.prop the proportion of zeros in the rows of the matrix.
#' @return A list with the following components: \item{R}{The sampled
#' correlation matrix} \item{L}{The loading matrix} 
#' @examples
#' R <- rcormat(m = 10)$R
#' eigen(R)
#' @export rcormat
rcormat <- function(m, r = 3L, sparse.prop = 0.5){
  fnums <- rnorm(m*r)
  L <- matrix(fnums, ncol = r, nrow = m)
  nzero <- sample(x = 0:(r-1), size = m, prob = c(sparse.prop, rep((1 - sparse.prop)/(r-1), 
              r-1)), replace = TRUE)
  for(i in 1L:m){
    if(nzero[i]) L[i,sample(1:r, nzero[i])] <- 0
  }
  L <- sweep(L, MARGIN = 1, STATS = sqrt(apply(L^2, 1, sum)), FUN = "/")  
  R <- tcrossprod(L)
  return(list(R = R, L = L))
}
