##################################################################
####
#### Find RIP-ordering of cliques of chordal (triangulated)
#### undirected graph
####
#### Based on Algorithm 4.11 in Steffen et all (the yellow book)
####
#### Known issues: Should check that amat specifies TUG;
#### possibly controlled by forceCheck argument
####
##################################################################

rip <- function(object, root=NULL, nLevels=NULL){
  UseMethod("rip")
}

rip.default <- function(object, root=NULL, nLevels=NULL){
    cls <- match.arg(class( object ),
                     c("graphNEL","igraph","matrix","dgCMatrix"))
    switch(cls,
           "graphNEL" ={
               if (graph::edgemode(object)=="directed"){
                   stop("Graph must be undirected")
               }
               ripMAT(graphNEL2dgCMatrix(object), root=root, nLevels=nLevels)
           },
           "igraph"   ={
               if (igraph::is.directed(object)){
                   stop("Graph must be undirected")
               }
               ripMAT(igraph::get.adjacency(object), root=root, nLevels=nLevels)
           },
           "dgCMatrix"=,
           "matrix"   ={
               ripMAT(object, root=root, nLevels=nLevels)
           })
}

ripMAT <- function(amat, root=NULL, nLevels=rep(2, ncol(amat))){

  ## mcs.idx: The enumeration of nodes in vn
  ## so 1,2,8 ... means that nodes 1,2,8... come first in the elimination
  ## mcs.vn: corresponding ordering of nodes

  mcs.vn <- mcsMAT(amat, root=root)
  if (length(mcs.vn)==0)
    return( NULL )

  cq.vn   <- maxCliqueMAT(amat)[[1]]
  ncq     <- length(cq.vn)
  vn      <- colnames(amat)

  if (ncq>1){
    ## cq.idx       <- lapplyV2I(cq.vn, vn)
    cq.mcs.idx   <- lapplyV2I(cq.vn, mcs.vn)
    max.no       <- sapply(cq.mcs.idx, function(cc) max(cc))
    ## max.no: The highest number assigned to a clique by MCS
    cq.mcs.idx   <- cq.mcs.idx[ order(max.no) ]
    ## cq.mcs.idx: Order the cliques by the largest number
    cq.mcs.vn    <- lapplyI2V( cq.mcs.idx, mcs.vn )

    sp.idx      <- vector("list", ncq)
    sp.idx[[1]] <- integer(0)
    pa.idx      <- rep.int(0L, ncq)

    ## Use: cq.mcs.idx
    for (ii in 2:ncq){
      paset  <- unlist(cq.mcs.idx[ 1:(ii-1L) ])
      isect  <- intersectPrim( cq.mcs.idx[[ii]], paset )
      sp.idx[[ii]] <- isect
      if (length(isect)){
        for (kk in (ii-1):1){
          if (subsetof( isect, cq.mcs.idx[[kk]]) ){
            pa.idx[ii] <- kk
            break()
          }
        }
      }
    }
    nodes <- mcs.vn
    sp.vn <- lapplyI2V(sp.idx, mcs.vn)
    child <- match(seq_along(cq.mcs.idx), pa.idx)
    host  <- rep.int(0L, length(nodes))
    ll    <- lapply(cq.mcs.vn, match, nodes)
    for (ii in seq_along(ll)){
      host[ll[[ii]]] <- ii
    }
  } else { ## The graph has only one clique!
    cq.mcs.vn <- list(vn)
    nodes  <- vn
    pa.idx <- 0L
    sp.vn  <- list(character(0))
    child  <- NA
    host  <- rep.int(1L, length(nodes))
  }

  .getChildList <- function(parents){
    vv  <- 1:length(parents)
    ft  <- unname(cbind(parents, vv))
    ft  <- ft[ft[, 1] != 0, , drop = FALSE]
    vch <- lapply(vv, function(ff) ft[ft[,1]==ff,2])
    names(vch) <- vv
    vch
  }
  ## FIXME: Create ftM2vch and ftM2vpar from the above

  rip3 <-
    structure(list(nodes       = nodes,
                   cliques     = cq.mcs.vn,
                   separators  = sp.vn,
                   parents     = pa.idx,
                   children    = child,
                   host        = host,
                   nLevels     = nLevels,
                   #createGraph = .createJTreeGraph,
                   childList   = .getChildList(pa.idx)
                   ),
              class="ripOrder")

  return(rip3)
}

print.ripOrder <- function(x, ...){
  idx <- 1:length(x$cliques)
  cat("cliques\n")
  mapply(function(xx,ii) cat(" ",ii,":",paste(xx, collapse=' '),"\n"), x$cliques, idx)

  cat("separators\n")
  mapply(function(xx,ii) cat(" ",ii,":",paste(xx, collapse=' '),"\n"), x$separators, idx)

  cat("parents\n")
  mapply(function(xx,ii) cat(" ",ii,":",paste(xx, collapse=' '),"\n"), x$pa, idx)

#  cat("Children\n")
#  mapply(function(xx,ii) cat(" ",ii,paste(xx, collapse=' '),"\n"), x$ch, idx)
}


plot.ripOrder <- function(x,...){

    .createJTreeGraph <- function(rip){
        if (length(rip$cliques)>1){
            ft <-cbind(rip$parents, 1:length(rip$parents))
            ft <- ft[ft[,1]!=0,, drop=FALSE]
            V <- seq_along(rip$parents)
            if (nrow(ft)==0){
                jt <- new("graphNEL", nodes = as.character(V), edgemode = "undirected")
            } else {
                jt <- graph::ftM2graphNEL(ft, V=as.character(V), edgemode="undirected")
            }
        } else {
            jt <- new("graphNEL", nodes = "1", edgemode = "undirected")
        }
        return(jt)
    }

    plot( .createJTreeGraph( x ) )
    ##plot(x$createGraph(x))
}

.rip2dag<-function (rip) {
  if (length(rip$cliques) > 1) {
    ft <- cbind(rip$parents, 1:length(rip$parents))
    ft <- ft[ft[, 1] != 0, , drop = FALSE]
    V  <- seq_along(rip$parents)
    if (nrow(ft) == 0) {
      jt <- new("graphNEL", nodes = as.character(V), edgemode = "directed")
    } else {
      jt <- graph::ftM2graphNEL(ft, V = as.character(V), edgemode = "directed")
    }
  } else {
    jt <- new("graphNEL", nodes = "1", edgemode = "directed")
  }
  return(jt)
}





