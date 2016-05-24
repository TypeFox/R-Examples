#' Converts to json representation
#' 
#' This is a nice function that takes a graph structure
#'  and converts it to a json format for use on the web.
#' @param graph An object of type \code{igraph}
#' @param file A file to write the output to (default missing)
#' @return A textual json representation of the graph
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
to_json <- function( graph, file ) {
  
  if( !inherits(graph,"popgraph"))
    stop("This function requires a popgraph object to function")
  
  # make vector 
  quotify <- function( df ){
    keys <- names(df)
    ret <- "["
    cols <- ncol(df)
    for(i in 1:nrow(df) ) {
      row <- "{"
      for( j in 1:ncol(df)){
        row <- paste(row,"\"",keys[j],"\":",sep="")
        valsep <- ifelse( is.numeric( df[i,j]), "", "\"" )
        row <- paste(row,valsep,df[i,j],valsep,sep="")
        if( j < ncol(df) )
          row <- paste( row, ",", sep="")
      }
      ret <- paste(ret, row,"}",sep="")
      
      if( i < nrow(df))
        ret <- paste(ret,", ",sep="")
      
    }
    ret <- paste( ret, ']',sep="")
    return(ret)
  }
  
  
  # do the nodes
  node.attr.names <- list.vertex.attributes( graph )
  if( !("name" %in% node.attr.names) )
    stop("Vertices are indexed by the property 'name' and your graph does not have one...")
  nodes <- data.frame( name=rep("libby",length(V(graph))) )
  for( attr in node.attr.names )
    nodes[[attr]] <- get.vertex.attribute( graph, attr )
  if( !("group" %in% names(nodes) ) )
    nodes$group <- "All"
  nodestr <- quotify(nodes)
  
  # make the edges
  K <- length(E(graph))
  edgedf <- data.frame( source=rep(1,K), target=rep(1,K) )
  if( "weight" %in% list.edge.attributes(graph))
    graph <- set.edge.attribute(graph,"weight", value=5)
  wts <- as.matrix(get.adjacency(graph, attr="weight"))
  idx <- 1
  N <- length(V(graph))
  for( i in 1:N){
    for( j in i:N) {
      if(wts[i,j] > 0 ) {
        edgedf$source[idx] <- (i-1)
        edgedf$target[idx] <- (j-1)
        idx <- idx + 1
      }
    }
  }
  edgestr <- quotify(edgedf)
  
  ret <- paste("var myjson = '{ \"nodes\":", nodestr, ", \"links\":",edgestr,"}';", sep="")
  
  if( !missing(file) )
    write(ret,file=file)
  else
    return( ret )
}