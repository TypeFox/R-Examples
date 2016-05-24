#' Returns distance congruence between the two graphs
#' 
#' This function makes the shortest path matrices for both
#'  graphs and determines the correlation between pairwise
#'  distance.
#' @param graph1 An object of type \code{igraph} or \code{popgraph}
#' @param graph2 An object of type \code{igraph} or \code{popgraph}
#' @param method An option on how congruence is to be estimated
#'  possible values are 'distance' (a measure of similarity in 
#'  separation of nodes independent of connectivity, the default) 
#'  'structural' a measure of similarity in actual edges, and 
#'  'combinatorial' a combinatorial measure of similarity.
#' @return A non-parametric rank sum test
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
test_congruence <- function( graph1, graph2, method=c("distance","combinatorial")[1] ) {
  cong.nodes <- intersect( V(graph1)$name , V(graph2)$name )
  
  if( is.null(cong.nodes) )
    stop("There appear to be no nodes in common between these two graphs")
  
  A <- induced.subgraph(graph1, vids=cong.nodes )
  B <- induced.subgraph(graph2, vids=cong.nodes )
  
  # do the distance congruence
  if( method == "distance" ) {
    Adis <- shortest.paths( A )
    Bdis <- shortest.paths( B )
    Adis <- Adis[cong.nodes, cong.nodes]
    Bdis <- Bdis[cong.nodes, cong.nodes]
    
    distances.graph1 <- Adis[lower.tri( Adis )]
    distances.graph2 <- Bdis[lower.tri( Bdis )]
    fit <- cor.test( distances.graph1, distances.graph2 )
  }
  
  
  # look at the structural congruence
  else if( method=="structural") {
    
    Aa <- as.matrix( get.adjacency( A ) )
    Ab <- as.matrix( get.adjacency( B ) )
    
    a11 <- 0.5 * (sum(Aa==1 & Ab==1) )
    a22 <- 0.5 * (sum(Aa==0 & Ab==0) - length(diag(Aa) ) )
    a12 <- 0.5 * (sum(Aa==0 & Ab==1) )
    a21 <- 0.5 * (sum(Aa==1 & Ab==0) )
    stop("Not Implemented Yet")
  }
  
  
  # combinatorial congruence
  else if( method=="combinatorial"){
    mA <- length( E(graph1) )
    mB <- length( E(graph2) )
    mC <- length( E(congruence_topology( graph1, graph2 ) ) ) 
    N <- length( V(graph1) )
    mMax <- N*(N-1)/2
    ell = mA + mB - mMax
    i <- max(0,ell):min(mA,mB)
    
    p.top <- choose( mA,mC ) * choose( mMax-mA, mB-mC )
    p.bot <- sum( choose(mA,i) * choose( mMax-mA, mB-i ))
    fit <- sum( p.top/p.bot )
    names(fit) <- "CDF"
  }
  
  # throw nothing back
  else
    fit <- NULL

  
  return( fit )
}