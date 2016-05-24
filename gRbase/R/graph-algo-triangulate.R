
triangulateMAT <- function(amat, nLevels=rep(2, ncol(amat)), ...){
    if (is.null(nLevels))
        nLevels=rep( 2, ncol(amat) )
    triangulateMAT_( amat, nLevels )
}


## FIXME: triangulate: Need clever choice of matrix-representation
## FIXME: (Sparse/dense)
triangulate.default <- function(object, nLevels=NULL, result=NULL, check=TRUE, ...)
{
    cls <- match.arg(class( object ),
                     c("graphNEL","igraph","matrix","dgCMatrix"))
    if (is.null( result ))
        result <- cls

    mm <- coerceGraph( object, "matrix" )
    if ( !is.UGMAT(mm) )
        stop("Graph must be undirected\n")

    if (check){
        if ( length(mcsMAT(mm)) > 0)
            coerceGraph(mm, result)
        else {
            mm <- triangulateMAT( mm, nLevels=nLevels )
            coerceGraph(mm, result)
        }
    } else {
        mm <- triangulateMAT( mm, nLevels=nLevels )
        coerceGraph(mm, result)
    }
}


















