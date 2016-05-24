#' Test modularity hypothesis
#'
#' Tests modularity hypothesis using cor.matrix matrix and trait groupings
#' @param cor.matrix Correlation matrix
#' @param modularity.hipot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hipot[i,j] == 1, trait i is in module j.
#' @param permutations Number of permutations, to be passed to MantelCor
#' @param MHI Indicates if test should use Modularity Hypothesis Index instead of AVG Ratio
#' @return Returns mantel correlation and associated probability for each modularity hypothesis, along with AVG+, AVG-, AVG Ratio for each module.
#' A total hypothesis combining all hypotesis is also tested.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MantelCor}}
#' @export
#' @rdname TestModularity
#' @references Porto, Arthur, Felipe B. Oliveira, Leila T. Shirai, Valderes Conto, and Gabriel Marroig. 2009. "The Evolution of Modularity in the Mammalian Skull I: Morphological Integration Patterns and Magnitudes." Evolutionary Biology 36 (1): 118-35. doi:10.1007/s11692-008-9038-3.
#' @examples
#' cor.matrix <- RandomMatrix(10)
#' rand.hipots <- matrix(sample(c(1, 0), 30, replace=TRUE), 10, 3)
#' mod.test <- TestModularity(cor.matrix, rand.hipots)
#' 
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#' cov.mod.test <- TestModularity(cov.matrix, rand.hipots, MHI = TRUE)
#' nosize.cov.mod.test <- TestModularity(RemoveSize(cov.matrix), rand.hipots, MHI = TRUE)
#' @keywords mantel
#' @keywords modularity
TestModularity <- function (cor.matrix, modularity.hipot, permutations = 100, MHI = FALSE) {
  m.hip.list <- CreateHipotMatrix(as.matrix(modularity.hipot))
  if(is.null(colnames(modularity.hipot))) colnames(modularity.hipot) <- 1:dim (modularity.hipot) [2]
  names(m.hip.list) <- c(colnames (modularity.hipot),"Full Integration")
  output <- MantelCor (m.hip.list, cor.matrix, permutations = permutations, mod = TRUE, MHI = MHI)
  names(output)[1] <- 'hypothesis'
  return (output)
}

#' @export
#' @rdname TestModularity
CreateHipotMatrix <- function(modularity.hipot) {
  num.hip <- dim (modularity.hipot) [2]
  num.traits <- dim (modularity.hipot) [1]
  m.hip.list <- alply(modularity.hipot, 2, function(x) outer(x, x))
  m.hip.list[[num.hip+1]] <- matrix(as.integer (as.logical (Reduce ("+", m.hip.list[1:num.hip]))),
                                    num.traits, num.traits, byrow=T)
  return(m.hip.list[1:(num.hip+1)])
}
