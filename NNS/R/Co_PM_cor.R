#' Co-Partial Moments Higher Dimension Correlation
#'
#' Determines higher dimension correlation coefficients based on degree 0 co-partial moments.
#' @param A data.frame of variables.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Beyond Correlation: Using the Elements of Variance for Conditional Means and Probabilities"  \url{http://ssrn.com/abstract=2745308}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100); z<-rnorm(100)
#' A<-data.frame(x,y,z)
#' \dontrun{Co.PM.cor(A)}
#' @export


Co.PM.cor <- function (A){

  n=ncol(A)

  A_upm = t(apply(A, 1, function(x) x>colMeans(A)))
  A_lpm = t(apply(A, 1, function(x) x<colMeans(A)))


  CO_upm = sum(apply(A_upm, 1, prod))/(length(A_upm)/n)
  CO_lpm = sum(apply(A_lpm, 1, prod))/(length(A_lpm)/n)


  observed = CO_upm+CO_lpm
  independence = 2*(.5^n)


  return((observed-independence)/(1-independence)) }
