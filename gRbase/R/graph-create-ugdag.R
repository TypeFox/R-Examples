##################################################################
####
#### Create undirected graphs or DAGs from graph specification
####
##################################################################

####################
## Undirected graphs
####################

ug <- function(...,result="NEL"){
  ugList(list(...), result=result)
}

ugList <- function(x, result="NEL"){
    result <- match.arg(result, c("matrix","Matrix","dgCMatrix","igraph","NEL","graphNEL"))
    x   <- unlist(lapply(x, function(g) rhsf2list(g)), recursive=FALSE)
    vn  <- unique.default( unlist(x, use.names=FALSE) )

    switch(result,
           "NEL"      =,
           "graphNEL" ={ugList2graphNEL(x, vn)},
           "Matrix"   =,
           "dgCMatrix"={ugList2dgCMatrix(x, vn)},
           "matrix"   ={ugList2matrix(x, vn)},
           "igraph"   ={
               gg <- igraph::igraph.from.graphNEL(ugList2graphNEL(x, vn))
               igraph::V(gg)$label <- igraph::V(gg)$name
               gg
           })
}


###########################
## Directed acyclic graphs
###########################

dag <- function(...,result="NEL", forceCheck=FALSE){
  dagList(list(...), result=result, forceCheck=forceCheck)
}

## dagList: forceCheck not implemented
dagList <- function(x, result="NEL", forceCheck=FALSE){
    result <- match.arg(result, c("matrix","Matrix","dgCMatrix","igraph","NEL","graphNEL"))
    x   <- unlist(lapply(x, function(g) rhsf2list(g)), recursive=FALSE)
    vn  <- unique(unlist(x))

    out <- switch(result,
                  "NEL"       =,
                  "graphNEL"  = {dagList2graphNEL(x, vn)},
                  "Matrix"    =,
                  "dgCMatrix" = {dagList2dgCMatrix(x, vn)},
                  "matrix"    = {dagList2matrix(x, vn)},
                  "igraph"    = {
                      gg <- igraph::igraph.from.graphNEL(dagList2graphNEL(x, vn))
                      igraph::V(gg)$label <- igraph::V(gg)$name
                      gg
                  })
    if (forceCheck){
        if( length( topoSort( out )) == 0){
            stop("In dag/dagList: Graph is not a DAG", call.=FALSE)
        }
    }
    out
}



