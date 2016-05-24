## ###############################################################
##
## Get list of edges in graph and list of edges not in graph
##
## FIXME: Check that it works on undirected and DAGs
##
## ###############################################################

edgeList <- function(object, matrix=FALSE)
  UseMethod("edgeList")

edgeList.default <- function(object, matrix=FALSE){
    cls <- match.arg(class( object ),
                     c("graphNEL","matrix","dgCMatrix"))
    switch(cls,
           "graphNEL"={edgeListMAT(graphNEL2M(object), matrix=matrix)},
           "dgCMatrix"=,
           "matrix"={edgeListMAT( object, matrix=matrix )})
}

edgeListMAT <- function(adjmat, matrix=FALSE){
    ans <-
        if (issymMAT_(adjmat)){
            symMAT2ftM_( adjmat )
        } else {
            MAT2ftM_( adjmat )
    }

    di <- dim(ans)
    ans <- colnames(adjmat)[ans]
    dim(ans) <- di

    if (!matrix){
        rowmat2list(ans)
    } else {
        ans
    }
}


nonEdgeList <- function(object, matrix=FALSE)
  UseMethod("nonEdgeList")

nonEdgeList.default <- function(object, matrix=FALSE){
    cls <- match.arg(class( object ),
                     c("graphNEL","matrix","dgCMatrix"))
    switch(cls,
           "graphNEL"={nonEdgeListMAT(graphNEL2M(object), matrix=matrix)},
           "dgCMatrix"=,
           "matrix"={nonEdgeListMAT( object, matrix=matrix )})

}

nonEdgeListMAT <- function(adjmat, matrix=FALSE){
    if (!issymMAT_( adjmat )) stop("'adjmat' must be symmetric")
    if (class(adjmat) == "dgCMatrix"){
        adjmat <- as( ((-1*adjmat) + 1), "dgCMatrix")
    } else {
        adjmat <- -1 * adjmat + 1
    }
    edgeListMAT( adjmat, matrix=matrix)
}


##
## vpar implemented for graphNEL, matrix and Matrix
##

vpar <- function(object, getv=TRUE, forceCheck=TRUE){
  UseMethod("vpar")
}

vparMAT <- function(object, getv=TRUE, forceCheck=TRUE){
  if (forceCheck && !is.adjMAT(object))
    stop("Matrix is not adjacency matrix... \n")
  if(forceCheck && isSymmetric(object))
    stop("Graph is undirected; (v,pa(v)) does not exist...\n")

  vn <- rownames(object)
  ##idx <- seq_along(vn)
  ans <- vector("list", length(vn))
  if (getv){
      for (jj in seq_along(vn)) {
          ans[[jj]] <- vn[c(jj, which(object[, jj]!=0))]
      }
  } else {
      for (jj in seq_along(vn)) {
          ans[[jj]] <- vn[object[, jj] != 0]
      }
  }
  names(ans) <- vn
  ans
}

vpar.graphNEL <- function(object, getv=TRUE, forceCheck=TRUE){
    if (forceCheck && graph::edgemode(object)=="undirected")
        stop("Graph is undirected; (v,pa(v)) does not exist...\n")

    ee <- graph::edges(object)
    vn <- names(ee)
    tf <- do.call(rbind, # matrix in to-from form
                  lapply(1:length(ee),
                         function(ii) names2pairs( ee[[ii]], vn[ii],
                                                  sort=FALSE, result="matrix")))

    ans <- lapply(1:length(vn), function(ii) c(vn[ii], tf[tf[,1]==vn[ii],2]))
    names(ans) <- vn
    if (!getv)
        ans<-lapply(ans, function(x)x[-1])
    return(ans)
}

vpar.Matrix <- vpar.matrix <- vparMAT


## #############################################################
##
## Get the cliques of an undirected graph
##
## #############################################################

## FIXME: Should check that it is undirected.
maxCliqueMAT <- function(amat){
  vn <- dimnames(amat)[[2L]]
  em <- t.default( MAT2ftM_( amat ) )
  maxClique(nodes=vn, edgeMat=em)
}

## FIXME: getCliques.graphNEL; graphNEL2dgCMatrix
## FIXME: -> should be graphNEL2adjMAT combined with an
## FIXME: -> intelligent choice of representation

getCliques <- function(object){
    UseMethod("getCliques")
}

getCliques.graphNEL <- function(object){
    maxCliqueMAT( graphNEL2dgCMatrix(object) )[[1]]
}

getCliques.default <- function(object){
    maxCliqueMAT(object)[[1]]
}



##
## Generate a random dag
##
random_dag <- function(V, maxpar=3, wgt=0.1){
    V <- as.character(V)
    vparList <- vector("list", length(V))
    names(vparList) <- V
    for (ii in 1:length(V)){
        rest <- V[-(1:ii)]
        zz <- 0:(min(maxpar, length(rest))-1)
        if (min(zz)<0)
            zz <- 0
        pp <- wgt^zz
        npar <- sample(zz, 1, prob=pp)
        vparList[[ii]] <- c(V[ii], sample(rest, npar, replace=FALSE))
    }

    dg <- dagList(vparList)
    dg
}


##
## SHD version of DED's dual rep; based on faster set operations
##
dual.rep <- function(glist, S, minimal=TRUE) {
    ## S is total varset - often but by no means always given by unique(unlist(g.list))
    list.save <- list()
    ##if (length(glist)==0) list.save <- list(S)
    if (length(glist)==1 & is.logical(glist[[1]]))
        list.save <- list(S)
    else {
        for (v in 1:length(glist)) {
            m1 <- list.save
            if (minimal)
                m2 <- as.list( setdiffPrim(S, glist[[v]]) )
            else
                m2 <- as.list( glist[[v]] )

            if (v==1)
                list.save <- m2
            else {
                ##aaa <- unlist(lapply(m1, function(g)
                ##                     lapply(m2, union, g)),recursive=FALSE)
                aaa <- unlist(lapply(m1, function(g)
                                     lapply(m2, function(o){unique.default(c(o, g))})),
                              recursive=FALSE)
                list.save <- removeRedundant(aaa, FALSE)
            }
        }
        if (!minimal)
            list.save <- lapply(list.save, function(g) setdiffPrim(S, g))}
    list.save
}





