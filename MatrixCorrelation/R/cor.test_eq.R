#' @title Test for no correlation between paired sampes
#'
#' @description Permutation test for squared Pearson correlation between to vectors of samples.
#'
#' @param x first \code{vector} to be compared (or two column \code{matrix/data.frame}).
#' @param y second \code{vector} to be compared (ommit if included in \code{x}).
#' @param B integer number of permutations, default = 10000.
#'
#' @details This is a convenience function combining \code{SMI} and \code{significant} for the
#' special case of vector vs vector comparisons. The nullhypothesis is that the correlation
#' between the vectors is +/-1, while significance signifies a deviance toward 0.
#'
#' @return A value indicating if the two input vectors are signficantly different.
#'
#' @author Kristian Hovde Liland
#'
#' @references Similarity of Matrices Index - Ulf Geir Indahl, Tormod Næs, Kristian Hovde Liland
#'
#' @seealso \code{\link{plot.SMI}} (print.SMI/summary.SMI), \code{\link{RV}} (RV2/RVadj), \code{\link{r1}} (r2/r3/r4/GCD), 
#' \code{\link{allCorrelations}} (matrix correlation comparison), \code{\link{PCAcv} (cross-validated PCA)}.
#'
#' @examples
#' a <- (1:5) + rnorm(5)
#' b <- (1:5) + rnorm(5)
#' cor.test_eq(a,b)
#'
#' @export
cor.test_eq <- function(x, y, B = 10000){
  # Split data if supplied as a matrix/data.frame
  if(missing(y)){
    y <- x[,2]
    x <- x[,1]
  }
  
  # Calculate SMI and significance
  smi <- SMI(x,y, 1,1)
  sig <- c(significant(smi, B = 10000))
  attr(sig, 'cor') <- cor(x,y)
  sig
}