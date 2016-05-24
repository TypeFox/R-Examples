
cnEdgeDistancePearson <- function(data, pert) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("'data' should be a matrix or data frame of cats")

  if(is.null(pert)) {
    warning("Perturbations are essential for estimating the pairwise causality")
  }

  if(is.data.frame(data)) {
    data <- as.matrix(t(data))
    if(!is.null(pert)) { 
      if(!is.data.frame(pert))
        stop("Perturbations should be a data frame")
      pert <- as.matrix(t(pert))
    }
  }
  
  r <- .categorizeSample(data, pert, object=NULL, ask=FALSE)
  data <- r$data
  pert <- r$pert
  numnodes <- dim(data)[1]
  numsamples <- dim(data)[2]
  nodenames <- rownames(data)
  
  mat <- .Call("ccnPearsonPairwise", 
                  data, pert, 
                  PACKAGE="sdnet")
  klmat <- matrix(mat, numnodes, numnodes)
  rownames(klmat)<-nodenames
  colnames(klmat)<-nodenames
  return(klmat)
}

.nodeChisq <- function(idroot, ppars, pcatlist, idx, problist, freqlist) {

  if(is.null(ppars) || length(idx) < 1) {

    if(length(problist) != length(pcatlist[[idroot]]) || length(pcatlist[[idroot]]) < 2) {
       return(list(chisq=0, df=0, nn=0))
    }
    
    nn <- sum(freqlist)
    if(is.numeric(nn) && nn > 0) {
      problist <- problist*nn
      diff <- problist - freqlist
      len <- length(pcatlist[[idroot]])
      chisq <- 0
      df <- 0
      for(i in 1:len) {
        if(problist[i]<=0)
          next
        chisq <- chisq + diff[i]*diff[i]/problist[i]
        df <- df + 1
      }
      ##cat(chisq, " ", df , "\n")
      return(list(chisq=chisq, df=df-1, nn=nn))
    }
    else
      return(list(chisq=0, df=0, nn=0))
  }
  idnode <- ppars[idx[1]]
  poutlist <- lapply(seq(1,length(pcatlist[[idnode]])),
                     function(cat) {
                       .nodeChisq(idroot, ppars, pcatlist, idx[-1], problist[[cat]], freqlist[[cat]])
                     })
  chisq <- 0
  df <- 0
  nn <- 0
  pcatlist[[idnode]]
  for(tt in poutlist) {
    chisq <- chisq + tt$chisq
    df <- df + tt$df
    nn <- nn + tt$nn
  }
  
  return(list(chisq=chisq, df=df, nn=nn))
}

setMethod("cnPearsonTest", signature("catNetwork"), 
          function(object, data) {

  if(!is(object, "catNetwork"))
    stop("Object should be catNetwork.")
 
  if(!is.matrix(data) && !is.data.frame(data))
    stop("'data' should be a matrix or data frame of cats")
  
  r <- .categorizeSample(data, NULL, object)
  data <- r$data
  
  if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
    stop("The number of nodes in  the object and data should be equal")
  
  rownames <- rownames(data)
  if(length(rownames) != object@numnodes)
    stop("The data rows should be named after the nodes of the object.")
  
  if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
    norder <- order(rownames)
    data <- data[norder,]
    rownames <- rownames(data)
    norder <- order(object@nodes)
    object <- cnReorderNodes(object, norder)
  }
  
  if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
    stop("The row names should correspond to the object nodes.")

  for(nnode in (1:object@numnodes)) {
    if(length(r$cats[[nnode]]) > length(object@cats[[nnode]])) {
       warning("Data has more cats than object: ",length(r$cats[[nnode]]), ", ", length(object@cats[[nnode]]), "\n")
       ## this can happen if the sample size is very small
       ## prune the sample
       ps <- data[nnode,]
       clen <- length(object@cats[[nnode]])
       for(i in 1:length(data[nnode,]))
         if(data[nnode,i]>clen) 
           data[nnode, i] <- as.integer(1+rbinom(1, clen-1, 0.5))
     }
  }
    
  numnodes <- dim(data)[1]
  numsamples <- dim(data)[2]

  chisq <- rep(0, object@numnodes)
  df <- rep(0, object@numnodes)
  nn <- rep(0, object@numnodes)
  
  for(nnode in (1:object@numnodes)) {
    pslot <- initSampleProb(nnode, 
                            object@pars[[nnode]],
                            object@cats,
                            seq(1,length(object@pars[[nnode]])))
    
    for(j in (1:numsamples)) {
      ps <- data[,j]
      pslot <- updateSampleProb(nnode, object@pars[[nnode]], object@cats,
                                seq(1,length(object@pars[[nnode]])), pslot, ps)
    }

    tt <- .nodeChisq(nnode, object@pars[[nnode]], object@cats, 
               seq(1,length(object@pars[[nnode]])),
               object@probs[[nnode]], pslot)
    chisq[nnode] <- tt$chisq
    df[nnode] <- tt$df
    nn[nnode] <- tt$nn
  }
  
  return(list(chisq=chisq, df = df, nn=nn))
} )

