## ##################################################################
##
## GRAPH PROPERTIES
##
## ##################################################################


isGraphical <- function(x){
    UseMethod("isGraphical")
}

isGraphical.default <- function( x ){
    if ((cls <- class( x )) %in% c("formula","list")){
        switch(cls,
               "formula"={
                   if (length(x)==3)
                       x <- formula(delete.response(terms(x)))
                   .isGraphical_glist( rhsf2list( x ) )
               },
               "list"={
                   if (!all(unlist(lapply(x, is.atomic))))
                       stop("'x' must be a list of atomic vectors")
                   if (length( x ) == 0)
                       stop("'x' must have positive length")
                   .isGraphical_glist( x )
               },
               stop("'x' must be a formula or a list of atomic vectors\n")
               )

    }
}

.isGraphical_glist <- function(x){
    vn <- unique( unlist(x) )
    amat <- glist2adjMAT(x, vn=vn)
    cliq <- maxCliqueMAT(amat)[[1]]
    all(unlist(lapply(cliq, function(sss) isin(x, sss))))
}

isDecomposable <- function(x){
    UseMethod("isDecomposable")
}

isDecomposable.default <- function( x ){
    if ((cls <- class( x )) %in% c("formula","list")){
        switch(cls,
               "formula"={
                   if (length(x)==3)
                       x <- formula(delete.response(terms(x)))
                   .isDecomposable_glist( rhsf2list( x ) )
               },
               "list"={
                   if (!all(unlist(lapply(x, is.atomic))))
                       stop("'x' must be a list of atomic vectors")
                   if (length( x ) == 0)
                       stop("'x' must have positive length")
                   .isDecomposable_glist( x )
               },
               stop("'x' must be a formula or a list of atomic vectors\n")
               )

    }
}

.isDecomposable_glist <- function(x){
    vn <- unique( unlist(x) )
    amat <- glist2adjMAT(x, vn=vn)
    cliq <- maxCliqueMAT(amat)[[1]]
    isg <- all(unlist(lapply(cliq, function(sss) isin(x, sss))))
    if (isg){
        length( mcsMAT( amat ) ) > 0
    } else
        FALSE
}

## Is model defined by <glist> graphical and strongly decomposable
## if discrete=NULL, then the check is just if the graph is decomposable
## Issues: Fails on the "empty graph".

isGSD_glist <- function(glist, vn=unique(unlist(glist)), discrete=NULL)
{
  amat <- glist2adjMAT(glist,vn=vn)
  cliq <- maxCliqueMAT(amat)[[1]]
  isg  <- all(unlist(lapply(cliq, function(sss) isin(glist, sss))))
  if (!isg){
    return(c(isg=FALSE, issd=FALSE))
  } else {
    return(c(isg=TRUE, issd=length(mcsmarkedMAT(amat,discrete=discrete)) > 0))
  }
}

properties_glist <- function(glist,
                             vn=unique(unlist(glist)),
                             amat=glist2adjMAT(glist,vn=vn),
                             cliq=maxCliqueMAT(amat)[[1]],discrete=NULL){

  isg <- all(unlist(lapply(cliq, function(sss) isin(glist, sss))))
  if (!isg){
    return(c(isg=FALSE, issd=FALSE))
  } else {
    return(c(isg=TRUE, issd=length(mcsmarkedMAT(amat,discrete=discrete)) > 0))
  }
}


##
## is.DAG, is.UG, is.TUG implemented for graphNEL, matrix and Matrix
##

is.adjMAT <- function(x){
    if (class(x) %in% c("matrix","dgCMatrix"))
        isadjMAT_( x )
    else
        FALSE
}

## ######################################

is.DAG <- function(object){UseMethod("is.DAG")}

is.DAG.graphNEL <- function(object){
    is.DAGMAT( graphNEL2dgCMatrix(object) )
}

is.DAG.default <- function( object ){
    if (class(object) %in% c("matrix","dgCMatrix"))
        isdagMAT_( object )
    else
        stop("'object' must be a matrix")
}

is.DAGMAT <- function(object){
    isdagMAT_( object )
}

#is.DAG.matrix <- is.DAG.Matrix <- is.DAGMAT

## ######################################

is.UG <- function(object){
    UseMethod("is.UG")
}

is.UG.graphNEL <- function(object){
    isugMAT_( graphNEL2MAT( object ) )
}

is.UG.default <- function( object ){
    if (class(object) %in% c("matrix","dgCMatrix"))
        isugMAT_( object )
    else
        stop("'object' must be a matrix")
}

is.UGMAT <- function(object){
    isugMAT_(object)
}

## ######################################

is.TUG <- function(object){
  UseMethod("is.TUG")
}

is.TUG.graphNEL <- function(object){
    z <- graphNEL2MAT( object )
    if (!isugMAT_( z ))
        FALSE
    else
        length( ripMAT( z ) )>0
}

is.TUG.default <- function(object){
    if (class(object) %in% c("matrix","dgCMatrix")){
        if (isugMAT_( object ))
            length(mcsMAT(object))>0
        else
            FALSE
    } else
        stop("'object' must be a matrix")
}

is.TUGMAT <- function(object){
    isugMAT_(object) && length(mcsMAT(object))>0
}

## ######################################

is.DG <- function(object){
    UseMethod("is.DG")
}

is.DG.graphNEL <- function(object){
    is.DGMAT( graphNEL2MAT(object) )
}

is.DG.default <- function(object){
    if (class(object) %in% c("matrix","dgCMatrix")){
        if (isadjMAT_(object))
            sum(object * t(object))==0
        else
            FALSE
    } else
        FALSE
}

is.DGMAT <- function(object){
    if (!is.adjMAT(object)) stop("Matrix is not adjacency matrix...\n")
    sum(object * t(object))==0
}



















## is.DG.matrix <- is.DG.Matrix <- is.DGMAT

## is.UG.graphNEL <- function(object){
##     is.UGMAT( graphNEL2dgCMatrix(object) )
## }

## is.UGMAT <- function(object){
##     isugMAT_( object )
## }

## is.UG.matrix <- is.UG.Matrix <- is.UGMAT

## .is.DAGMAT <- function(object){
##     if (!is.adjMAT(object)) stop("Matrix is not adjacency matrix...\n")
##     length(topoSort(object))>0
## }

## .is.UGMAT <- function(object){
##     if (!is.adjMAT(object)) stop("Matrix is not adjacency matrix...\n")
##     isSymmetric(object)
## }


## isGraphical.formula <- function(x){
##     if (length(x)==3)
##         x <- formula(delete.response(terms(x)))
##     .isGraphical_glist( rhsf2list( x ) )
## }

## isGraphical.list <- function(x){
##     if (!all(unlist(lapply(x, is.atomic))))
##         stop("'x' must be a list of atomic vectors")
##     if (length( x ) == 0)
##         stop("'x' must have positive length")
##     .isGraphical_glist( x )
## }

## isDecomposable.formula <- function(x){
##     if (length(x)==3)
##         x <- formula(delete.response(terms(x)))
##     .isDecomposable_glist( rhsf2list( x ) )
## }

## isDecomposable.list <- function(x){
##     if (!all(unlist(lapply(x, is.atomic))))
##         stop("'x' must be a list of atomic vectors")
##     if (length( x ) == 0)
##         stop("'x' must have positive length")
##     .isDecomposable_glist( x )
## }


## .is.adjMAT <- function(x){
##   res <- FALSE
##   if (inherits(x, c("matrix","Matrix"))){
##     d <- dim(x)
##     if (d[1L]==d[2L]){
##       v <- 1:d[1L]
##       if( all(x[cbind(v, v)]==0) ){
##         res <- TRUE
##       }
##     }
##   }
##   res
## }
