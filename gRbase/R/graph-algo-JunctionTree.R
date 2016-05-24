jTree <- function(object, ...){
  UseMethod("jTree")
}

jTree.default <-  function(object, nLevels = NULL, ...){
    cls <- match.arg(class( object ),
                     c("graphNEL","igraph","matrix","dgCMatrix"))

    switch(cls,
           "graphNEL" ={jTreeMAT( graphNEL2dgCMatrix( object ), nLevels=nLevels, ... )},
           "igraph"   ={jTreeMAT( igraph::get.adjacency( object ), nLevels=nLevels, ... )},
           "dgCMatrix"=,
           "matrix"   ={jTreeMAT( object, nLevels=nLevels, ... )})
}


jTreeMAT <- function(amat, nLevels=rep(2,ncol(amat)), ...){
  tug  <- triangulateMAT( amat, nLevels=nLevels, result="Matrix", ... )
  ripMAT( tug, nLevels=nLevels )
}


## Copying
junctionTree          <- jTree
junctionTree.default  <- jTree.default
junctionTreeMAT       <- jTreeMAT


## jTree.graphNEL <- function(object, nLevels = NULL, ...){
##   jTreeMAT( graphNEL2dgCMatrix( object ), nLevels=nLevels, ... )
## }

## jTree.igraph <- function(object, nLevels = NULL, ...){
##   jTreeMAT( igraph::get.adjacency( object ), nLevels=nLevels, ... )
## }

## jTree.matrix <- function(object, nLevels = NULL, ...){
##   jTreeMAT( object, nLevels=nLevels, ... )
## }

## jTree.Matrix <- jTree.matrix

## junctionTree.graphNEL <- jTree.graphNEL
## junctionTree.igraph   <- jTree.igraph
## junctionTree.matrix   <- jTree.matrix
## junctionTree.Matrix   <- jTree.Matrix






























## jTree.graphNEL <- function(object, method  = "mcwh",
##                   nLevels = rep(2,length(nodes(object))), ...){

##   method <- match.arg(tolower(method),c("mcwh","r"))
##   tug        <- triangulate(object, method=method, nLevels=nLevels, result="Matrix")
##   val        <- ripMAT(tug,nLevels=nLevels)
##   return(val)
## }

## jTree.matrix <- function(object, method  = "mcwh",
##                          nLevels = rep(2,ncol(object)), ...){

##   method <- match.arg(tolower(method),c("mcwh","r"))
##   tug        <- triangulateMAT(object, method=method, nLevels=nLevels, result="Matrix")
##   val        <- ripMAT(tug,nLevels=nLevels)
##   return(val)
## }

## jTree.Matrix <- function(object, method  = "mcwh",
##                          nLevels = rep(2,ncol(object)), ...){

##   method <- match.arg(tolower(method),c("mcwh","r"))
##   tug        <- triangulateMAT(object, method=method, nLevels=nLevels, result="Matrix")
##   val        <- ripMAT(tug,nLevels=nLevels)
##   return(val)
## }

## jTree.igraph <- function(object,  method  = "mcwh",
##                   nLevels = rep(2,length(V(object))), ...){

##   method <- match.arg(tolower(method),c("mcwh","r"))
##   tug        <- triangulate(object, method=method, nLevels=nLevels, result="matrix")
##   val        <- ripMAT(tug,nLevels=nLevels)
##   return(val)
## }





