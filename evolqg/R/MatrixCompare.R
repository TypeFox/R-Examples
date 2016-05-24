#' Matrix Compare
#'
#' Compare two matrices using all available methods. Currently RandomSkewers, MantelCor, KrzCor and PCASimilarity
#'
#' @param cov.x covariance or correlation matrix
#' @param cov.y covariance or correlation matrix
#' @param id name of the comparison column
#' @return data.frame of comparisons
#' @export
#' @examples
#'
#' cov.x = RandomMatrix(10, 1, 1, 10)
#' cov.y = RandomMatrix(10, 1, 10, 20)
#' MatrixCompare(cov.x, cov.y)

MatrixCompare <- function(cov.x, cov.y, id = '.id'){
  RS <- RandomSkewers(cov.x, cov.y)[1:2]
  Mantel <- MantelCor(cov2cor(cov.x), cov2cor(cov.y))
  Krz <- c(KrzCor(cov.x, cov.y), NA)
  PCA <- c(PCAsimilarity(cov.x, cov.y), NA)
  output <- data.frame(rbind(RS, Mantel, Krz, PCA))
  output[id] = rownames(output)
  return(output)
}
