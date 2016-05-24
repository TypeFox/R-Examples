#' Prune a graph
#'
#' @param g sg object
#' @param level pruning level
#' @param verbose verbosity
#'
#' @export

prune.sg<-function(g, level=1, verbose=FALSE) {
  if(!is(g,"sg")) stop("g not sg object.")
  if(is.null(level)) return(g)
  if(level<=0)return(g)

  g <- sg2sym(g)

  edges <- prune_c(g$edges, level, verbose)

  as.sg(edges,type=g$type, pars=g$parameters,
        note=c(g$note, paste("pruned with level=", as.integer(level),sep="")))
}
