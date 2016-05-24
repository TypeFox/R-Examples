# In graph theory, a topological sort or topological ordering of a
# directed acyclic graph (DAG) is a linear ordering of its nodes in
# which each node comes before all nodes to which it has outbound edges.
# Every DAG has one or more topological sorts. If such ordering can
# not be found then the graph has cycles
#
# Input:  list of vectors of the form (v,pa(v))
# Output: vector with ordering
#
# should perhaps be called dagTopoSort


topoSort <- function(object, index=FALSE){
  UseMethod("topoSort")
}

topoSort.default <- function(object, index=FALSE){
    cls <- match.arg(class( object ),
                     c("graphNEL","igraph","matrix","dgCMatrix"))
    switch(cls,
           "graphNEL" ={topoSortMAT(graphNEL2dgCMatrix(object), index=index) },
           "igraph"   ={topoSortMAT(igraph::get.adjacency(object), index=index) },
           "dgCMatrix"=,
           "matrix"   ={topoSortMAT(object, index=index)} )
}

topoSortMAT <- function(XX_, index=FALSE){
    ans <- topoSortMAT_( XX_ )
    if (index){
        if (ans[1]!=-1){
            ans
        } else {
            -1L
        }
    } else {
        if (ans[1]!=-1){
            colnames(XX_)[ans]
        } else {
            character(0)
        }
    }
}

topoSort_vparList<- function(vpaL){
    topoSort(vpaList2adjMAT(vpaL, result="Matrix"))
}



## topoSort.graphNEL<- function(object, index=FALSE){
##   topoSortMAT(as(object,"Matrix"), index=index)
## }



## topoSort.matrix <- topoSort.Matrix <- function(object, index=FALSE){
##   topoSortMAT(object, index=index)
## }

## topoSortMAT <- function(XX_, index=FALSE){
##   if (inherits(XX_, "Matrix")){
##     ans <- .Call("gRbase_topoSortMAT_sp", XX_ ,package="gRbase")
##   } else {
##     if (inherits(XX_, "matrix")){
##       ans <- .Call("gRbase_topoSortMAT_st", XX_ ,package="gRbase")
##     } else {
##       stop("'XX_' must be a matrix or a sparse matrix (a 'dgCMatrix')")
##     }
##   }
##   if (index){
##     if (ans[1]!=-1){
##       ans
##     } else {
##       -1L
##     }
##   } else {
##     if (ans[1]!=-1){
##       colnames(XX_)[ans]
##     } else {
##       character(0)
##     }
##   }
## }
