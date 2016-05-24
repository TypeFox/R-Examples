#' Alpha Repetability
#'
#' Calculates the matrix repeatability using the equation in Cheverud 1996
#' Quantitative genetic analysis of cranial morphology in the cotton-top
#' (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
#'
#'
#' @param cor.matrix Correlation matrix
#' @param sample.size Sample size used in matrix estimation
#' @return Alpha repeatability for correlation matrix
#' @keywords repeatatability
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{BootstrapRep}}
#' @author Diogo Melo, Guilherme Garcia
#' @references Cheverud 1996 Quantitative genetic analysis of cranial morphology in the cotton-top
#'   (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
#' @export
#' @examples
#' #For single matrices
#' cor.matrix <- RandomMatrix(10)
#' AlphaRep(cor.matrix, 10)
#' AlphaRep(cor.matrix, 100)

#' #For many matrices
#' mat.list <- RandomMatrix(10, 100)
#' sample.sizes <- floor(runif(100, 20, 50))
#' unlist(Map(AlphaRep, mat.list, sample.sizes))

AlphaRep <- function (cor.matrix, sample.size){
  if(sum(diag(cor.matrix)) != dim(cor.matrix)[1])
    stop("Matrices do not appear to be correlation matrices.")
  vec <- cor.matrix[lower.tri(cor.matrix)]
  var.erro <- (1 - mean(vec)^2)/(sample.size-2)
  var.vec <- var(vec)
  return((var.vec - var.erro)/var.vec)
}
