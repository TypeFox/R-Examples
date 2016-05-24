#' Cut a quantitative variable in \eqn{n} equal parts
#'
#' Cuts a quantitative variable in \eqn{n} equal parts
#'
#' @param x a numeric vector
#' @param n numeric, the number of parts: 2 to cut according to the median, and so on...
#' @param \dots other arguments to be passed in \code{\link{cut}}
#' @return A factor vector
#' @author Hugo Varet
#' @examples
#' cut_quanti(cgd$height, 3)

cut_quanti=function(x,n,...){
  out=cut(x,breaks=quantile(x,probs=seq(0,1,length=n+1),na.rm=TRUE),include.lowest=TRUE,...)
  return(out)
}

#cut_quanti(rnorm(100),3,dig.lab=1)

