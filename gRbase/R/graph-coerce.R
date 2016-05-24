##############################################################
####
#### Coercion between graphNEL, igraph and matrix
####
##############################################################

setOldClass("igraph")

matrix2igraph <- function(x){
    if (isSymmetric(x)){
        gg <- igraph::graph.adjacency(x, mode="undirected")
    } else {
        gg <- igraph::graph.adjacency(x, mode="directed")
    }
    igraph::V(gg)$label <- igraph::V(gg)$name <- colnames(x)
    gg
}

graphNEL2igraph <- function(x){
    gg <- igraph::igraph.from.graphNEL(x)
    igraph::V(gg)$label <- igraph::V(gg)$name
    gg
}








## From graphNEL
## -------------
setAs("graphNEL", "igraph",      function(from) graphNEL2igraph(from))
setAs("graphNEL", "matrix",      function(from) graphNEL2M(from, result="matrix"))
setAs("graphNEL", "Matrix",      function(from) graphNEL2M(from, result="Matrix"))
setAs("graphNEL", "dgCMatrix",   function(from) graphNEL2M(from, result="Matrix"))

## From matrix
## -----------
setAs("matrix", "igraph", function(from) matrix2igraph(from))

## matrix -> graphNEL : is in graph package (I guess)
## matrix -> dgCMatrix: is in Matrix package. Should be used
## matrix -> Matrix : is in Matrix package but care should be taken
## because the output can be of different types

## From Matrix
## -----------
setAs("Matrix", "igraph", function(from){ matrix2igraph( as.matrix( from )) })

# Matrix -> graphNEL : in the graph package (I guess)
# Matrix -> matrix : in the Matrix package

## From igraph
## -----------
setAs("igraph",   "graphNEL",    function(from) igraph::igraph.to.graphNEL(from))
setAs("igraph",   "matrix",      function(from) as(igraph::get.adjacency(from),"matrix"))
setAs("igraph",   "Matrix",      function(from) MAT2dgCMatrix(igraph::get.adjacency(from)))
setAs("igraph",   "dgCMatrix",   function(from) MAT2dgCMatrix(igraph::get.adjacency(from)))



## #####################################################################
##
## The coerceGraph methods are mentioned in GMwR and therefore they
## must be kept alive.
##
## #####################################################################

coerceGraph <- function(object, result){
  UseMethod("coerceGraph")
}

coerceGraph.graphNEL <- function(object, result){
  result <- match.arg(result, c("graphNEL","matrix","dgCMatrix","Matrix","igraph"))
  switch(result,
         "graphNEL"={object},
         "igraph"  ={gg <- igraph::igraph.from.graphNEL(object)
                     igraph::V(gg)$label <- igraph::V(gg)$name
                     gg
                   },
         "matrix" =,
         "Matrix" =,
         "dgCMatrix"={
           graphNEL2M(object, result=result)
         }
         )
}

coerceGraph.matrix <- function(object, result){
  result <- match.arg(result, c("graphNEL","matrix","dgCMatrix","Matrix","igraph"))
  switch(result,
         "graphNEL" ={ as(object,"graphNEL")},
         "igraph"   ={ matrix2igraph(object)},
         "matrix"   ={ object },
         "Matrix"   =,
         "dgCMatrix"={ matrix2dgCMatrix( object )})
}

coerceGraph.dgCMatrix <- function(object, result){
  result <- match.arg(result, c("graphNEL","igraph","matrix","dgCMatrix","Matrix"))
  switch(result,
         "graphNEL" ={ as(object,"graphNEL")},
         "igraph"   ={ matrix2igraph(dgCMatrix2matrix(object))},
         "matrix"   ={ dgCMatrix2matrix( object )},
         "Matrix"   =,
         "dgCMatrix"={ object  })
}


coerceGraph.igraph <- function(object, result){
  result <- match.arg(result, c("graphNEL","matrix","dgCMatrix","Matrix","igraph"))
  switch(result,
         "graphNEL"={ igraph::igraph.to.graphNEL(object)},
         "igraph"  ={ object},
         "matrix"  ={ as(igraph::get.adjacency(object),"matrix")},
         "Matrix"  =,
         "dgCMatrix"={ MAT2dgCMatrix(igraph::get.adjacency(object))}
         )
}



### xxx2yyy

ugList2graphNEL<- function(gset, vn=NULL){
    if ( is.null(vn) )
        vn <- unique.default( unlist(gset, use.names=FALSE) )
    zzz <- lapply(gset, function(xx) names2pairs(xx, sort=TRUE, result="matrix"))
    ftM <- do.call(rbind, zzz)
    if ( nrow(ftM) > 0 ){
        tofrom <- unique(rowmat2list(ftM))
        fff <- do.call(rbind, tofrom)
        graph::ftM2graphNEL(fff, V=as.character(vn), edgemode="undirected")
    } else {
        new("graphNEL", nodes=as.character(vn), edgemode="undirected")
    }
}

dagList2graphNEL<- function(gset, vn=NULL){
    if ( is.null(vn) )
        vn <- unique.default( unlist(gset, use.names=FALSE) )
    zzz <- lapply(gset, function(xx) names2pairs(xx[1],xx[-1],
                                                    sort=FALSE, result="matrix"))
    ftM <- do.call(rbind, zzz)
    if (nrow(ftM)>0){
        tfL <- unique(rowmat2list(ftM))
        ftM <- do.call(rbind,tfL)[,2:1,drop=FALSE]
        graph::ftM2graphNEL(ftM, V=as.character(vn),
                            edgemode="directed")
    } else {
        new("graphNEL", nodes=as.character(vn), edgemode="directed")
    }
}




##################################################
##
## Convert between matrix and dgCMatrix
##
##################################################

MAT2matrix <- function(x){
    .check.that.input.is.matrix(x)
    switch( class(x),
           "matrix"    ={x},
           "dgCMatrix" ={dgCMatrix2matrix(x)})
}

MAT2dgCMatrix <- function(x){
    .check.that.input.is.matrix(x)
    switch( class(x),
           "matrix"    ={matrix2dgCMatrix(x)},
           "dgCMatrix" ={x})
}

##################################################
##
## Convert list of generators to adjacency matrix
##
##################################################

## glist: A list of vectors of the form (v, pa1, pa2, ... pan)
vpaList2adjMAT <- function(glist, vn=unique(unlist(glist)), result="matrix"){
    result <- match.arg(result, c("matrix", "Matrix", "dgCMatrix"))
    switch(result,
           "Matrix"=,
           "dgCMatrix" = {dagList2dgCMatrix( glist, vn )},
           "matrix"    = {dagList2matrix( glist, vn )}  )
}

## glist: A list of vectors of the form (v1, v2, ... vn)
glist2adjMAT <- function(glist, vn=unique(unlist(glist)), result="matrix"){
    result <- match.arg(result, c("matrix","Matrix","dgCMatrix"))
    switch(result,
           "Matrix"=,
           "dgCMatrix" = {ugList2dgCMatrix( glist, vn )},
           "matrix"    = {ugList2matrix( glist, vn )}  )
}

## adjList : named list as returned by graph::edges( )
adjList2adjMAT <- function(adjList, result="matrix"){
    result <- match.arg(result, c("matrix", "Matrix", "dgCMatrix"))
    switch(result,
           "matrix"   = {adjList2matrix( adjList )},
           "Matrix"   = ,
           "dgCMatrix"= {adjList2dgCMatrix( adjList )})
}

adjList2M <- function( x, result="matrix"){
    adjList2adjMAT(x, result=result)
}



##
## graphNEL 2 something
##

graphNEL2M <- function(object, result="matrix"){
    if( class(object) != "graphNEL" )
        stop("'object' must be a graphNEL object...")
    adjList2adjMAT( graph::edges(object), result=result )
}

## FIXME graphNEL2adjMAT used by HydeNet package; I do not use it.
graphNEL2adjMAT <- graphNEL2M
as.adjMAT       <- graphNEL2M

## Never used
graphNEL2matrix    <- function(object){
    graphNEL2M(object, result="matrix")
}
## Used a lot
graphNEL2dgCMatrix <- function(object){
    graphNEL2M(object, result="Matrix")
}

graphNEL2MAT <- function(object, limit=100){
    if( class(object) != "graphNEL" )
        stop("'object' must be a graphNEL object...")

    result <-
        if ( length( graph::nodes(object) ) > limit )
            "dgCMatrix" else "matrix"

    adjList2M( graph::edges(object), result=result )
}


## vpaL2tfM: (v,pa(v))-list 2 to-from-matrix
## FIXME vpaL2tfM: rename to vpaList2ftM; used in topoSort
vpaL2tfM <- function(vpaL){
 eMat  <- lapply(vpaL, function(xx) names2pairs(xx[1], xx[-1],
                                                 sort = FALSE, result = "matrix"))
 do.call(rbind, eMat)
}

graphNEL2ftM <- function(object){
    if( class(object) != "graphNEL" )
        stop("'object' must be a graphNEL object...")
    adjList2ftM(graph::edges(object))
}

graphNEL2tfM <- function(object){
    if( class(object) != "graphNEL" )
        stop("'object' must be a graphNEL object...")
    adjList2tfM(graph::edges(object))
}

## -----------
ugList2M <- function(x, result="matrix"){
    ## glist2adjMAT <- function(glist, vn=unique(unlist(glist)), result="matrix")
    result <- match.arg(result, c("matrix","Matrix","dgCMatrix"))
    vn <- unique.default(unlist(x), use.names=FALSE)
    switch(result,
           "Matrix"=,
           "dgCMatrix" = {ugList2dgCMatrix( x, vn )},
           "matrix"    = {ugList2matrix( x, vn )}  )
}

## -----------
dagList2M <- function(x, result="matrix"){
    ## vpaList2adjMAT(x, result=result)
    result <- match.arg(result, c("matrix", "Matrix", "dgCMatrix"))
    vn <- unique.default(unlist(x), use.names=FALSE)
    switch(result,
           "Matrix"=,
           "dgCMatrix" = {dagList2dgCMatrix( x, vn )},
           "matrix"    = {dagList2matrix( x, vn )}  )
}


##
## Matrix 2 something
##

M2adjList <- function(x){
    .check.that.input.is.matrix(x)
    vn <- colnames(x)
    if (!isadjMAT_(x))
        stop("'x' is not an adjacency matrix\n")
    r  <- rowmat2list(x)
    i  <- lapply(r, function(z) which(z!=0))
    out <- lapply(i, function(j) vn[j])
    names(out) <- vn
    out
}

M2ugList <- function(x){
    ## FIXME: M2ugList: Need a check for undirectedness
    .check.that.input.is.matrix(x)
    maxCliqueMAT(x)[[1]]
}

M2graphNEL <- function(x){
    .check.that.input.is.matrix(x)
    as(x, "graphNEL")
}

M2dagList <- function(x){
    .check.that.input.is.matrix(x)
    vn <- colnames(x)
    c  <- colmat2list(x)
    i  <- lapply(c, function(z) which(z!=0))
    i  <- lapply(1:length(vn), function(j) c(j, i[[j]]))
    out <- lapply(i, function(j) vn[j])
    ##names(out) <- vn
    out
}

.check.that.input.is.matrix <- function(x){
    if ( !(class(x)=="matrix" || class(x)=="dgCMatrix") )
        stop("Input must be a matrix or a dgCMatrix\n")
}


ug2dag <- function(object){

    if (class(object) != "graphNEL")
        stop("Object 'object' must be a graphNEL")
    if (graph::edgemode(object) != "undirected")
        stop("Graph must have undirected edges")
    if (length( m <- mcs(object) )==0)
        stop("Graph is not chordal")

    adjList  <- graph::adj(object, m)
    vparList <- vector("list", length(m))
    names(vparList) <- m

    vparList[[1]] <- m[1]
    if (length(m) > 1){
        for (i in 2:length(m)){
            vparList[[ i ]] <- c(m[ i ],
                                intersectPrim(adjList[[ i ]], m[ 1:i ]))
        }
    }

    dg <- dagList(vparList)
    dg
}


#' .eliminationOrder <- function(gg){
#'   is.acyc <- TRUE
#'   ### amat <- as.adjmat(gg)
#'   amat <- as.adjMAT(gg)
#'   elorder <- NULL

#'   repeat{
#'     idx <- which(rowSums(amat)==0)
#'     if (!length(idx)){
#'       return(NULL)
#'     }
#'     elorder <- c(elorder, idx)
#'     amat <- amat[-idx,-idx]

#'     if(all(c(0,0)==dim(amat))){
#'       break()
#'     }
#'   }
#'   names(rev(elorder))
#' }






















## Represent list of sets in a matrix...
## FIXME: glist2setMAT: Used in gRain 1.2-3, but not in gRain 1.2-4
## FIXME: should be deleted for next release
glist2setMAT <- function(glist,vn=unique(unlist(glist))){
  amat <- matrix(0, nrow=length(glist), ncol = length(vn))
  colnames(amat) <- vn
  for (i in 1:length(glist)){
    amat[i, glist[[i]] ] <- 1
  }
  amat
}




#' genL2M <- function( x, result="matrix"){
#'     ##glist2adjMAT <- function(glist, vn=unique(unlist(glist)), result="matrix")
#'     ugList2M(x, result=result)
#' }

#' vpaL2M <- function(x, result="matrix"){
#'     vpaList2adjMAT(x, result=result)
#' }




