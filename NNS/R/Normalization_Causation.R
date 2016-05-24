#' VN.norm.caus
#'
#' VN.norm function with reduced output specifically for VN.caus routine
#' @param A Matrix of variables.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' A<-cbind(x,y)
#' \dontrun{VN.norm.caus(A)}
#' @export


VN.norm.caus <- function(A) {
  m  <- colMeans(A)
  RG <- m %o% (1/m)
  scales <- colMeans(RG * abs(cor(A)))
  A_Normalized <- t(t(A) * scales)
  return(A_Normalized)}
