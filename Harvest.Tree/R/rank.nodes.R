#' Ranking of nodes
#' 
#' Rank harvested node by lower p value
#' @param harfunc.object an object of class "harfunc"
#' @return the ranked harvest nodes
#' @export



rank.nodes <- function(harfunc.object){
  har.rule <- harfunc.object$har.rule
  lowci.har<- unlist(lapply(har.rule, function(x) {Bi.test(x$total, x$active)}))
  lowci.n <- rank(-lowci.har)
  return(lowci.n)
}