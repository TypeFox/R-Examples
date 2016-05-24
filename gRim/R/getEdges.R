#########################################################################
###
### Get edges *in* or *not in* a graph
###
### Input: Various formats
### Output: A p * 2 matrix
###
### getInEdges: if type="decomposable" then the edges returned are those
### edges e for which the graph is decomposable if e is removed. If
### type="unrestricted" then all edges are returned
###
### getOutEdges: if type="decomposable" then the edges returned are those
### edges e for which the graph is decomposable if e is added. If
### type="unrestricted" then all edges are returned
###
### Known issues: None
###
##########################################################################

getEdges <- function(object,type="unrestricted", ingraph=TRUE, discrete=NULL, ...){
  UseMethod("getEdges")
}

getEdges.list <- function(object,type="unrestricted", ingraph=TRUE, discrete=NULL, ...){
  if (ingraph){
    getInEdgesMAT(glist2adjMAT(object), type=type, discrete=discrete)
  } else {
    getOutEdgesMAT(glist2adjMAT(object), type=type, discrete=discrete)
  }
}

getEdges.graphNEL <- function(object,type="unrestricted", ingraph=TRUE, discrete=NULL, ...){
  if (ingraph){
    getInEdgesMAT(as.adjMAT(object), type=type, discrete=discrete)
  } else {
    getOutEdgesMAT(as.adjMAT(object), type=type, discrete=discrete)
  }
}

getEdges.matrix <- function(object,type="unrestricted", ingraph=TRUE, discrete=NULL, ...){
  if (ingraph){
    getInEdgesMAT(object, type=type, discrete=discrete)
  } else {
    getOutEdgesMAT(object, type=type, discrete=discrete)
  }
}

##
## Convenience functions...
##
getInEdges <- function(object, type="unrestricted", discrete=NULL, ...){
  getEdges(object, type=type, ingraph=TRUE, discrete=discrete, ...)
}

getOutEdges <- function(object, type="unrestricted", discrete=NULL, ...){
  getEdges(object, type=type, ingraph=FALSE, discrete=discrete, ...)
}


##########################################################################

getInEdgesMAT <- function(adjmat, type="unrestricted", discrete=NULL, ...){
  type <- match.arg(type, c("unrestricted", "decomposable"))
  emat <- edgeListMAT(adjmat, matrix=TRUE) #;  print(emat)
  if (type=="decomposable"){
      idx  <- vector("logical", nrow(emat))    
      for (ii in seq_len(nrow(emat))){
        ed <- emat[ii,] #;        print(ed)
        adjmat[ed[1],ed[2]] <- adjmat[ed[2],ed[1]] <- 0
        idx[ii] <- length(mcsmarkedMAT(adjmat,discrete=discrete))>0
        adjmat[ed[1],ed[2]] <- adjmat[ed[2],ed[1]] <- 1
      }
      emat <- emat[idx,,drop=FALSE]
    }
  emat
}

getOutEdgesMAT <- function(adjmat, type="unrestricted", discrete=NULL, ...){
  type <- match.arg(type, c("unrestricted", "decomposable"))
  emat <- nonEdgeListMAT(adjmat, matrix=TRUE)
  if (type=="decomposable"){
      idx  <- vector("logical", nrow(emat))    
      for (ii in seq_len(nrow(emat))){
        ed <- emat[ii,]
        adjmat[ed[1],ed[2]] <- adjmat[ed[2],ed[1]] <- 1
        idx[ii] <- length(mcsmarkedMAT(adjmat, discrete=discrete))>0
        adjmat[ed[1],ed[2]] <- adjmat[ed[2],ed[1]] <- 0
      }
      emat <- emat[idx,,drop=FALSE]
    }
  emat
}
