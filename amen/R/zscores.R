#' rank-based z-scores
#' 
#' Computes the normal scores corresponding to the ranks of a data vector
#' 
#' 
#' @usage zscores(y)
#' @param y a numeric vector
#' @return a numeric vector
#' @author Peter Hoff
#' @export zscores
zscores<-function(y)
{
 qnorm( rank(y,na.last="keep",ties.method="average")/(1+sum(!is.na(y))) )
}

