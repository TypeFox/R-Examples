#' Converts graph to a pgraph file format
#' 
#' This is a simple function that takes the graph and
#'  converts it into a *.pgraph file for visualization 
#'  in other software.
#' @param graph An object of type \code{popgraph}
#' @param file The name and location of where the *.pgraph file is to be saved.
#'  If ommitted, this function will return a single text file.
#' @return Nothing if passed a file or the raw text of the *.pgraph file if
#'  you do not provide a file object.
#' @export
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
to_pgraph <- function( graph, file ) {
  if( !is(graph,"popgraph") )
    stop("This function only works using a popgraph object.")

  K <- length(V(graph))
  L <- length(E(graph))
  
  if( !("name" %in% list.vertex.attributes(graph)))
    V(graph)$name <- paste("Node",1:K,sep="-")
  if( !("size" %in% list.vertex.attributes(graph)))
    V(graph)$size <- 1
  if( !("color" %in% list.vertex.attributes(graph)))
    V(graph)$color <- 1
  
  if( !("weight" %in% list.edge.attributes(graph) ) )
    E(graph)$weight <- 1
  
  pgraph_text <- paste( K, "\t", L, sep="")
  
  for( i in 1:K )
    pgraph_text <- append( pgraph_text, paste( V(graph)$name[i], "\t", V(graph)$size[i], "\t", V(graph)$color[i], sep="") )
  
  A <- get.adjacency(graph,attr="weight")
  for( i in 1:K )
    for( j in i:K)
      if( A[i,j] > 0 )
        pgraph_text <- append( pgraph_text, paste( V(graph)$name[i],"\t",V(graph)$name[j],"\t",A[i,j],sep="") )
  
  pgraph_text <- paste( pgraph_text, collapse="\n")
  
  
                               
                               
  if( !missing(file) ) {
    write(pgraph_text,file)
    invisible(pgraph_text)
  }
  else
    return( pgraph_text )
}

