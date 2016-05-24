#' Converts \code{popgraph} to \code{data.frame} based upon node attributes
#' 
#' This is a quick conversion of vertex attributes to a \code{data.frame}, essentially
#'  the reverse operation as \code{decorate_graph} function.
#' @param x The \code{popgraph} to grab stuff from.
#' @param ... Ignored (generally).
#' @return An object of type \code{data.frame} with all the node attributes.
#' @author Rodney J. Dyer \email{rjdyer@@vcu.edu}
#' @export
to_data.frame <- function( x, ... ){
  if( !is(x,"popgraph") & !(is(x,"igraph")))
    stop("What are you passing to to_data.frame()?")
  
  cols <- list.vertex.attributes( x )
  ret <- data.frame( vertex.id=seq(1,length(V(x))))
  if( length(cols) > 0 ) {
    ret[[cols[1]]] <- get.vertex.attribute(x,name=cols[1])
    for( i in 2:length(x))
      ret[[cols[i]]] <- get.vertex.attribute(x, name=cols[i])  
  }
  if( !("vertex.id" %in% cols))
    ret$vertex.id <- NULL
  
  return(ret)
}