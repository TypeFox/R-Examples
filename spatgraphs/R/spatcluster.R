#' Compute the connected components of a graph
#'
#' @param x sg-object
#' @param verbose print info
#' @param sym force symmetry of edges
#'
#' @export

spatcluster <- function(x, verbose=TRUE, sym =FALSE){
  if(sym) x<-sg2sym(x)

  clusters<- spatcluster_c(x$edges, verbose)

  as.sgc(clusters, x$type, x$parameters)

}
