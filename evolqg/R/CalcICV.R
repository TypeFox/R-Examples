#' Calculates the ICV of a covariance matrix.
#'
#' Calculates the coefficient of variation of the eigenvalues of a covariance matrix, a measure of 
#' integration comparable to the R^2 in correlation matrices.
#'
#' @param cov.matrix Covariance matrix.
#' @return coefficient of variation of the eigenvalues of a covariance matrix
#' @export
#' @author Diogo Melo
#' @seealso \code{\link{CalcR2}}
#' @references Shirai, Leila T, and Gabriel Marroig. 2010. "Skull Modularity in Neotropical Marsupials and Monkeys: Size Variation and Evolutionary Constraint and Flexibility." Journal of Experimental Zoology Part B: Molecular and Developmental Evolution 314 B (8): 663-83. doi:10.1002/jez.b.21367.
#' @references Porto, Arthur, Leila Teruko Shirai, Felipe Bandoni de Oliveira, and Gabriel Marroig. 2013. "Size Variation, Growth Strategies, and the Evolution of Modularity in the Mammalian Skull." Evolution 67 (July): 3305-22. doi:10.1111/evo.12177.
#' @examples
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#' CalcICV(cov.matrix)
#' @keywords covariance
#' @keywords integration

CalcICV <- function (cov.matrix){
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(sum(diag(cov.matrix)) == dim(cov.matrix)[1]) warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used for ICV.")
  eVals <- eigen(cov.matrix, only.values = TRUE)$values
  ICV <- sd(eVals)/mean(eVals)
  return (ICV)
}