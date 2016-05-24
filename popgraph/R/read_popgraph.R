#' Reads in a population graph text file 
#' 
#' This function imports a popgraph from the older text format
#' @param file The path to the file that is to be saved.
#' @param sep The column separator in the file. By default it is a 
#'  tab but sometimes a space or other object may be used.
#' @return A fully created popgraph file (e.g., an igraph object
#'  with an extra class property)
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
read.popgraph <- function( file, sep="\t" ) { 
    
  # load in the raw stuff
  raw <- read.table( file,header=FALSE, stringsAsFactors=FALSE,sep=sep,fill=TRUE)
  
  
  # set up the adjacency matrix
  Knode <- as.numeric(raw[1,1])
  Kedge <- as.numeric(raw[1,2])
  
  # set up some colors
  names <- raw[2:(Knode+1),1]
  sizes <- as.numeric(raw[2:(Knode+1),2])
  colors <- rep("#FDAE61",Knode)
  
  A <- matrix( 0, ncol=Knode, nrow=Knode)
  
  i<-Knode+2
  j<-Knode+Kedge+1
  raw1 <- raw[i:j,]
  
  # go through the edges 
  for( i in 1:length(raw1[,1])){
    row <- raw1[i,]    
    fidx <- which( names == row[1,1] )
    tidx <- which( names == row[1,2] )
    wt <- as.numeric( row[1,3])
    if( length( fidx ) & length( tidx) ) 
      A[fidx,tidx] <- A[tidx,fidx] <- wt
  }
  
  rownames(A) <- colnames(A) <- names
  graph <- graph.adjacency( A, mode="undirected", weighted=TRUE)
  V(graph)$size <- sizes 
  V(graph)$color <- colors
  
  class(graph) <- c("igraph","popgraph")
  return( graph )
  
}