#' Robust and classical variation matrix
#' 
#' Estimates the variation matrix with robust methods.
#' 
#' The variation matrix is estimated for a given compositional data set.
#' Instead of using the classical standard deviations the
#' \code{\link[stats]{mad}} is used when parameter robust is set to TRUE.
#' 
#' @param x data frame or matrix with positive entries
#' @param robust if FALSE, standard measures are used.
#' @return The (robust) variation matrix.
#' @author Matthias Templ
#' @references Aitchison, J. (1986) \emph{The Statistical Analysis of
#' Compositional Data} Monographs on Statistics and Applied Probability.
#' Chapman \& Hall Ltd., London (UK). 416p.
#' @keywords multivariate robust
#' @export
#' @examples
#' 
#' data(expenditures)
#' variation(expenditures)
#' variation(expenditures, robust=FALSE)
#' 
`variation` <-
  function(x, robust=TRUE){
    rvars <- matrix(0, ncol=ncol(x), nrow=ncol(x))
    if(robust){
      for( i in 1:ncol(x)){
        for( j in 1:ncol(x)){
          if( i < j ){
            rvars[i,j] <- (mad(log(x[,i]/x[,j])))^2
            rvars[j,i] <- rvars[i,j]
          }
        }
      }
    } else{
      for( i in 1:ncol(x)){
        for( j in 1:ncol(x)){
          if( i < j ){ 
            rvars[i,j] <- (var(log(x[,i]/x[,j])))
            rvars[j,i] <- rvars[i,j]            
          }
        }
      }		
    }
    return(rvars) 
}
