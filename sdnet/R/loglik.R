
#########################################################################
# Categorical Network Class Methods
# Log-Likelihood

.nodeLikelihood <- function(idroot, ppars, pcatlist, idx, problist, psample) {
  if(is.null(ppars) || length(idx) < 1) {
    if(is.na(psample[idroot]))
      return(NA)
    ##cat(idroot," ", psample[idroot], "\n")
    if(psample[idroot] < 1 || psample[idroot] > length(pcatlist[[idroot]]))
      stop("Wrong sample\n")
    cat <- psample[idroot]
    ##cat(idroot, ": nodeLikelihood ", cat, ", ", problist[cat], "\n")
    return(problist[cat])
  }
  idnode <- ppars[idx[1]]
  if(idnode < 1 || idnode > length(psample))
    stop("Wrong sample")
  cat <- psample[idnode]
  if(is.na(cat))
    return(NA)
  ##cat(idroot," ", idnode, " cat = ", cat)
  return(.nodeLikelihood(idroot, ppars, pcatlist, idx[-1], problist[[cat]], psample))
}

.networkLikelihood <- function(object, data) {
  if(dim(data)[1] != object@numnodes)
    stop("Incompatible sample dimensions.\n")
  
  numSamples <- dim(data)[2]
  if(numSamples < 1)
    stop("No samples\n")

  r <- .categorizeSample(data, NULL, object)
  data <- r$data
  
  floglik <- 0
  for(j in (1:numSamples)) {
    ps <- data[,j]
    liklist <- sapply(seq(1,object@numnodes), function(nnode) {
      .nodeLikelihood(nnode, object@pars[[nnode]], object@cats, 
                      seq(1,length(object@pars[[nnode]])),
                      object@probs[[nnode]], ps)
    })
    liklist <- liklist[!is.na(liklist)]
    floglik <- floglik + mean(sapply(liklist, function(x) if(x > 0) log(x) else -Inf))
  }

  return(floglik)
}


setMethod("cnLoglik", c("catNetwork"), 
          function(object, data, pert=NULL, bysample=FALSE, softmode=FALSE) { 

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of cats")

            if(is.data.frame(data)) {
              data <- as.matrix(t(data))
              if(!is.null(pert)) { 
                if(!is.data.frame(pert))
                  stop("Perturbations should be a data frame")
                pert <- as.matrix(t(pert))
              }
            }

            if(softmode) {
              if(!is.numeric(data) || is.integer(data))
                stop("Numeric data is expected in soft mode")
              rownames <- object@nodes
              numnodes <- length(object@nodes)
              numSamples <- ncol(data)
            }
            else {

              if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
                stop("The number of nodes in  the object and data should be equal")
            
              rownames <- rownames(data)
              if(length(rownames) != object@numnodes)
                stop("The data rows/cols should be named after the nodes of the object.")
            
              if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
                norder <- order(rownames)
                ## keep the matrix format !!
                data <- as.matrix(data[norder,])
                if(!is.null(pert))
                  pert <- as.matrix(pert[norder,])
                rownames <- rownames(data)
                norder <- order(object@nodes)
                object <- cnReorderNodes(object, norder)
              }
              
              if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
                stop("The data names should correspond to the object nodes.")
              
              r <- .categorizeSample(data, NULL, object, ask=FALSE)
              data <- r$data
              if(object@maxcats < r$maxcats)
                stop("Data has more cats than the object")
              
              numnodes <- dim(data)[1]
              numSamples <- dim(data)[2]
            }

            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The data names should correspond to the object nodes.")
            
            loglik <- .Call("ccnLoglik", 
                            object, data, pert, bysample, 
                            PACKAGE="sdnet")
            if(bysample)
              return(loglik)
            return(sum(loglik))
          })


.nodeSampleLoglik <- function(nnode, parentSet, data, cats, maxcats) {

  numnodes <- dim(data)[1]
  numSamples <- dim(data)[2]

  nodeprob <- initSampleProb(nnode, parentSet, cats, seq(1,length(parentSet)))
  for(j in (1:numSamples)) {
    ps <- data[,j]
    ## increment frequency      
    nodeprob <- updateSampleProb(nnode, parentSet, cats,
                                 seq(1,length(parentSet)), nodeprob, ps)
  }

  nodeprob <- normalizeProbSlot(nnode, parentSet, cats,
                   seq(1,length(parentSet)), nodeprob)

  # calculate loglik (a kind of redundancy)
  floglik <- 0
  if(numSamples < 1)
    return(floglik)
  floglik <- NULL
  for(j in (1:numSamples)) {
    ps <- data[,j]
    lik <- .nodeLikelihood(nnode, parentSet, cats,
                      seq(1,length(parentSet)),
                      nodeprob, ps)
    if(is.na(lik))
       next
    ##cat(nnode, ": ", lik, "\n")
    floglik <- c(floglik, lik)
  }
  return(mean(log(floglik)))
}

cnNodeSampleLoglik <- function(node, pars, data, pert = NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame of node cats")

  if(is.numeric(data) && !is.integer(data))
      stop("only categorical data is accepted")
  
  if(is.matrix(data)) {
    numnodes <- dim(data)[1]
    numSamples <- dim(data)[2]
    nodenames <- rownames(data)
  }
  else {
    numnodes <- dim(data)[2]
    numSamples <- dim(data)[1]
    nodenames <- colnames(data)
  }
  if(!is.character(nodenames) || length(nodenames) != numnodes)
    stop("Named nodes are expected")
 
  if(numSamples < 1)
    stop("No samples are given\n")

  if(length(pars) == 0)
    pars <- NULL
   
  r <- .categorizeSample(data, pert)
  data <- r$data
  pert <- r$pert
  cats <- r$cats
  maxcats <- r$maxcats

  if(is.character(node))
    node <- sapply(node, function(nn) return(which(nodenames == nn)))
  if(!is.numeric(node) || sum(node < 1) || sum(node > numnodes))
    stop("Incorect node ",node)
  node <- as.integer(node)
  
  if(length(node)==1 && !is.list(pars)) {
    if(node < 1 || node > numnodes)
      stop("Incorrect node ", node)
    if(length(pars)>0 && is.character(pars))
      pars <- sapply(pars, function(pp) which(nodenames == pp))
    pars <- as.integer(pars)
    
    loglik <- -Inf
    if(!is.null(pert)) {      
      subdata <- data[, pert[node,]==0]
      numsubsamples <- nrow(subdata)
      if(numsubsamples > 0) {
        loglik <- .nodeSampleLoglik(node, pars, subdata, cats, maxcats)
      }
    }
    else
      loglik <- .nodeSampleLoglik(node, pars, data, cats, maxcats)
    return(loglik)
  }

  if(length(node) > 1 && is.list(pars)) {
    if(length(node) != length(pars))
      stop("`node' and `pars' lists should be of equal length")
    loglik <- rep(-Inf, length(node))
    for(i in 1:length(node)) {
      nod <- node[[i]]
      par <- pars[[i]]
      if(length(par)>0 && is.character(par))
        par <- sapply(par, function(pp) which(nodenames == pp))
      par <- as.integer(par)
      
      if(nod < 1 || nod > numnodes)
        stop("Incorrect node ", nod)
      if(!is.null(pert)) {
        subdata <- data[, pert[nod,]==0]
        numsubsamples <- nrow(subdata)
        if(numsubsamples > 0) {
          loglik[i] <- .nodeSampleLoglik(nod, par, subdata, cats, maxcats)
        }
      }
      else {
        loglik[i] <- .nodeSampleLoglik(nod, par, data, cats, maxcats)
      }
    }
    return(loglik)
  }

  stop("`node' should be either an integer or list of integers")
  return(0)
}

cnNodeSampleProb <- function(node, pars, data, pert = NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame of node cats")

  if(is.numeric(data) && !is.integer(data))
      stop("only categorical data is accepted")
  
  if(is.matrix(data)) {
    numnodes <- dim(data)[1]
    numSamples <- dim(data)[2]
    nodenames <- rownames(data)
  }
  else {
    numnodes <- dim(data)[2]
    numSamples <- dim(data)[1]
    nodenames <- colnames(data)
  }
  
  if(numSamples < 1)
    stop("No samples are given\n")
  if(length(node) != 1)
    stop("Only one node is expected\n")
  
  if(length(pars) == 0)
    pars <- NULL
  
  if(length(nodenames) > 0) {
    if(is.character(node))
      node <- which(nodenames == node)
    if(is.character(pars))
      pars <- sapply(pars, function(par) which(nodenames == par))
  }

  if(is.character(node))
    node <- sapply(node, function(nn) return(which(nodenames == nn)))
  node <- as.integer(node)
  if(!is.integer(node))
    stop("Not valid node")
  if(node < 1 || node > numnodes)
    stop("Incorrect node ", node)
  if(!is.null(pars)) {
    if(length(pars)>0 && is.character(pars))
      pars <- sapply(pars, function(pp) which(nodenames == pp))
    pars <- as.integer(pars)
    if(!is.integer(pars))
      stop("Not valid pars")
  }

  r <- .categorizeSample(data, pert)
  data <- r$data
  pert <- r$pert
  cats <- r$cats
  maxcats <- r$maxcats
  
  nodeprob <- initSampleProb(node, pars, cats, seq(1,length(pars)))
  for(j in (1:numSamples)) {
    ps <- data[,j]
    ## increment frequency
    nodeprob <- updateSampleProb(node, pars, cats,
                                 seq(1,length(pars)), nodeprob, ps)
  }
  nodeprob <- normalizeProbSlot(node, pars, cats,
                                seq(1,length(pars)), nodeprob)
  return(nodeprob)
}

setMethod("cnNodeLoglik", c("catNetwork"), 
          function(object, node, data, pert=NULL, softmode=FALSE, klmode=0) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of cats")

            if(is.data.frame(data)) {
              data <- as.matrix(t(data))
              if(!is.null(pert)) { 
                if(!is.data.frame(pert))
                  stop("Perturbations should be a data frame")
                pert <- as.matrix(t(pert))
              }
            }

            if(softmode) {
              if(!is.numeric(data) || is.integer(data))
                stop("Numeric data is expected in soft mode")
              rownames <- object@nodes
              numnodes <- length(object@nodes)
              numSamples <- ncol(data)
            }
            else {
              if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
                stop("The number of nodes in  the object and data should be equal")
            
              rownames <- rownames(data)
              if(length(rownames) != object@numnodes)
                stop("The data rows/cols should be named after the nodes of the object.")
              
              if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
                norder <- order(rownames)
                ## keep the matrix format !!
                data <- as.matrix(data[norder,])
                if(!is.null(pert))
                  pert <- as.matrix(pert[norder,])
                rownames <- rownames(data)
                norder <- order(object@nodes)
                object <- cnReorderNodes(object, norder)
                if(is.numeric(node))
                  node <- sapply(node, function(nn) return(which(norder == nn)))
              }

              r <- .categorizeSample(data, pert, object, nodeCats=NULL, ask=FALSE)
              data <- r$data
              pert <- r$pert
              if(object@maxcats < r$maxcats)
                stop("Data has more cats than the object")
              
              numnodes <- dim(data)[1]
              numSamples <- dim(data)[2]
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The data names should correspond to the object nodes.")

            if(is.character(node))
              node <- sapply(node, function(nn) return(which(rownames == nn)))
            if(!is.numeric(node) || sum(node < 1) || sum(node > object@numnodes))
              stop("Wrong node")
            node <- as.integer(node)

            if(TRUE) {
              loglik <- .Call("ccnNodeLoglik", 
                              object, node, data, pert, klmode, 
                              PACKAGE="sdnet")
              return(loglik)
            }
            else
              return(.nodeLikelihood(node,
                                     object@pars[[node]], object@cats, 
                                     seq(1,length(object@pars[[node]])),
                                     object@probs[[node]], data))
          })
