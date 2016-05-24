moralize <- function(object,...){
  UseMethod("moralize")
}
moralize.default <- function(object, result=NULL, ...)
{
    cls <- match.arg(class( object ),
                     c("graphNEL", "matrix", "dgCMatrix", "igraph"))
    if (is.null( result ))
        result <- cls

    switch(cls,
           "graphNEL" ={
               m <- graphNEL2M( object )
               if ( !is.DAGMAT( m ) )
                   stop("Graph must be directed")
               m <- moralizeMAT( m )
           },
           "dgCMatrix"=,
           "matrix"   ={
               if ( !is.DAGMAT( object ) )
                   stop("Graph must be directed")
               m <- moralizeMAT(object)
           },
           "igraph"   ={
               if (!igraph::is.directed( object ))
                   stop("Graph must be directed")
               if (is.null(igraph::V(object)$name))
                   igraph::V(object)$name <- igraph::V(object)
               m <- moralizeMAT(igraph::get.adjacency(object))
           })
    coerceGraph( m, result)
}




