#' Mean Squared Correlations
#'
#' Calculates the mean squared correlation of a covariance or correlation matrix. Measures integration.
#'
#' @param c.matrix Covariance or correlation matrix.
#' @return Mean squared value of off diagonal elements of correlation matrix
#' @export
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{Flexibility}}
#' @references Porto, Arthur, Felipe B. de Oliveira, Leila T. Shirai, Valderes de Conto, and Gabriel Marroig. 2009. "The Evolution of Modularity in the Mammalian Skull I: Morphological Integration Patterns and Magnitudes." Evolutionary Biology 36 (1): 118-35. doi:10.1007/s11692-008-9038-3.
#' @references Porto, Arthur, Leila Teruko Shirai, Felipe Bandoni de Oliveira, and Gabriel Marroig. 2013. "Size Variation, Growth Strategies, and the Evolution of Modularity in the Mammalian Skull." Evolution 67 (July): 3305-22. doi:10.1111/evo.12177.
#' @examples
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#' # both of the following calls are equivalent, 
#' # CalcR2() converts covariance matrices to correlation matrices internally
#' CalcR2(cov.matrix)
#' CalcR2(cov2cor(cov.matrix))
#' @keywords correlation
#' @keywords integration

CalcR2 <- function (c.matrix){
    cor.matrix = cov2cor(c.matrix)
    return (mean (cor.matrix [lower.tri (cor.matrix)]^2))
}
