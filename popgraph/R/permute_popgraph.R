#' Bootstraps indiviudals to see stability of graph topology.
#'
#' This function uses a permutation test to look at edge stability.  What we 
#'  do is resample individuals and re-estimate the topology several times. This
#'  provides an estimate of edge stability.
#' @param data The raw multivariate data as submitted to \code{popgraph}
#' @param groups The grouping of the data into nodes as submitted to \code{popgraph}
#' @param nboot The nubmer of times to bootstrap the individuals per group (default=50)
#' @param ... Other arguments to be passed to \code{popgraph}
#' @return A weighted graph where edge weights represent the proportion of times the 
#'  edge was found in the perumuted data sets.
#' @export
#' 
permute_popgraph <- function( data, groups, nboot=50, ...){
  if( !is(data,"matrix"))
    stop("Cannot use non-matrix data to make a graph, let alone bootstrap it...")
  if( nrow(data) != length(groups))
    stop("You need to have data of the same size to use this function.")
  if( !is(groups,"factor"))
    groups <- factor(groups)
  
  strata <- levels(groups)
  K <- length(strata)
  A <- matrix(0,K,K)
  df <- data.frame(data)
  rownames(A) <- colnames(A) <- levels(groups)
  df$Stratum <- as.numeric(as.factor(groups))
  nc <- ncol(df)-1
  sz <- as.numeric( table(df$Stratum))
  for( rep in 1:nboot){
    ndata <- as.matrix(df[strata( df, stratanames = "Stratum",size=sz,method="srswr")$ID_unit,1:nc]  )
    graph <- popgraph(ndata,groups)
    B <- to_matrix(graph,mode = "adjacency")
    B[ B!=0 ] <- 1
    A <- A+B
    cat(".")
    if( !(rep %% 50))
      cat(" [",rep,"/",nboot,"]\n",sep="")
  }
  A <- A/nboot
  orig <- popgraph( data, groups )
  origA <- to_matrix(orig,mode = "adjacency")
  A <- A * origA
  
  graph <- decorate_graph(orig, A )
  return(graph)
}