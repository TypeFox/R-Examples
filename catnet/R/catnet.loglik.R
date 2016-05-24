
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
  
  numsamples <- dim(data)[2]
  if(numsamples < 1)
    stop("No samples\n")

  r <- .categorizeSample(data, NULL, object)
  data <- r$data
  
  floglik <- 0
  for(j in (1:numsamples)) {
    ps <- data[,j]
    liklist <- sapply(seq(1,object@numnodes), function(nnode) {
      .nodeLikelihood(nnode, object@parents[[nnode]], object@categories, 
                      seq(1,length(object@parents[[nnode]])),
                      object@probabilities[[nnode]], ps)
    })
    liklist <- liklist[!is.na(liklist)]
    floglik <- floglik + mean(sapply(liklist, function(x) if(x > 0) log(x) else -Inf))
  }

  return(floglik)
}


setMethod("cnLoglik", c("catNetwork"), 
          function(object, data, perturbations=NULL, bysample=FALSE) { 

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of categories")

            if(is.data.frame(data)) {
              data <- as.matrix(t(data))
              if(!is.null(perturbations)) { 
                if(!is.data.frame(perturbations))
                  stop("Perturbations should be a data frame")
                perturbations <- as.matrix(t(perturbations))
              }
            }

            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
              stop("The number of nodes in  the object and data should be equal")
            
            rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows/cols should be named after the nodes of the object.")
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              ## keep the matrix format !!
              data <- as.matrix(data[norder,])
              if(!is.null(perturbations))
                perturbations <- as.matrix(perturbations[norder,])
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The data names should correspond to the object nodes.")

            r <- .categorizeSample(data, NULL, object, ask=FALSE)
            data <- r$data
            if(object@maxCategories < r$maxCategories)
              stop("Data has more categories than the object")
            
            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]
            
            loglik <- .Call("ccnLoglik", 
                            object, data, perturbations, bysample, 
                            PACKAGE="catnet")
            ##return(.networkLikelihood(object, data))
            ##loglik <- loglik[!is.nan(loglik)]
            if(bysample)
              return(loglik)
            return(sum(loglik))
          })


.nodeSampleLoglik <- function(nnode, parentSet, data, categories, maxCategories) {

  numnodes <- dim(data)[1]
  numsamples <- dim(data)[2]

  nodeprob <- initSampleProb(nnode, parentSet, categories, seq(1,length(parentSet)))
  for(j in (1:numsamples)) {
    ps <- data[,j]
    ## increment frequency      
    nodeprob <- updateSampleProb(nnode, parentSet, categories,
                                 seq(1,length(parentSet)), nodeprob, ps)
  }

  nodeprob <- normalizeProbSlot(nnode, parentSet, categories,
                   seq(1,length(parentSet)), nodeprob)

  # calculate likelihood (a kind of redundancy)
  floglik <- 0
  if(numsamples < 1)
    return(floglik)
  floglik <- NULL
  for(j in (1:numsamples)) {
    ps <- data[,j]
    lik <- .nodeLikelihood(nnode, parentSet, categories,
                      seq(1,length(parentSet)),
                      nodeprob, ps)
    if(is.na(lik))
       next
    ##cat(nnode, ": ", lik, "\n")
    floglik <- c(floglik, lik)
  }
  return(mean(log(floglik)))
}

cnNodeSampleLoglik <- function(node, parents, data, perturbations = NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame of node categories")

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

  if(length(parents) == 0)
    parents <- NULL
   
  r <- .categorizeSample(data, perturbations)
  data <- r$data
  perturbations <- r$perturbations
  categories <- r$categories
  maxCategories <- r$maxCategories

  if(is.character(node))
    node <- sapply(node, function(nn) return(which(nodenames == nn)))
  if(!is.numeric(node) || sum(node < 1) || sum(node > numnodes))
    stop("Incorect node ",node)
  node <- as.integer(node)
  
  if(length(node)==1 && !is.list(parents)) {
    if(node < 1 || node > numnodes)
      stop("Incorrect node ", node)
    if(length(parents)>0 && is.character(parents))
      parents <- sapply(parents, function(pp) which(nodenames == pp))
    parents <- as.integer(parents)
    
    loglik <- -Inf
    if(!is.null(perturbations)) {      
      subdata <- data[, perturbations[node,]==0]
      numsubsamples <- nrow(subdata)
      if(numsubsamples > 0) {
        loglik <- .nodeSampleLoglik(node, parents, subdata, categories, maxCategories)
      }
    }
    else
      loglik <- .nodeSampleLoglik(node, parents, data, categories, maxCategories)
    return(loglik)
  }

  if(length(node) > 1 && is.list(parents)) {
    if(length(node) != length(parents))
      stop("`node' and `parents' lists should be of equal length")
    loglik <- rep(-Inf, length(node))
    for(i in 1:length(node)) {
      nod <- node[[i]]
      par <- parents[[i]]
      if(length(par)>0 && is.character(par))
        par <- sapply(par, function(pp) which(nodenames == pp))
      par <- as.integer(par)
      
      if(nod < 1 || nod > numnodes)
        stop("Incorrect node ", nod)
      if(!is.null(perturbations)) {
        subdata <- data[, perturbations[nod,]==0]
        numsubsamples <- nrow(subdata)
        if(numsubsamples > 0) {
          loglik[i] <- .nodeSampleLoglik(nod, par, subdata, categories, maxCategories)
        }
      }
      else {
        loglik[i] <- .nodeSampleLoglik(nod, par, data, categories, maxCategories)
      }
    }
    return(loglik)
  }

  stop("`node' should be either an integer or list of integers")
  return(0)
}

cnNodeSampleProb <- function(node, parents, data, perturbations = NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame of node categories")

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
  if(length(node) != 1 || length(parents) != 1)
    stop("Only one node is expected\n")
  
  if(length(parents) == 0)
    parents <- NULL
  
  if(length(nodenames) > 0) {
    if(is.character(node))
      node <- which(nodenames == node)
    if(is.character(parents))
      parents <- sapply(parents, function(par) which(nodenames == par))
  }

  if(is.character(node))
    node <- sapply(node, function(nn) return(which(nodenames == nn)))
  node <- as.integer(node)
  if(!is.integer(node))
    stop("Not valid node")
  if(node < 1 || node > numnodes)
    stop("Incorrect node ", node)
  if(!is.null(parents)) {
    if(length(parents)>0 && is.character(parents))
      parents <- sapply(parents, function(pp) which(nodenames == pp))
    parents <- as.integer(parents)
    if(!is.integer(parents))
      stop("Not valid parents")
  }

  r <- .categorizeSample(data, perturbations)
  data <- r$data
  perturbations <- r$perturbations
  categories <- r$categories
  maxCategories <- r$maxCategories
  
  nodeprob <- initSampleProb(node, parents, categories, seq(1,length(parents)))
  for(j in (1:numSamples)) {
    ps <- data[,j]
    ## increment frequency
    nodeprob <- updateSampleProb(node, parents, categories,
                                 seq(1,length(parents)), nodeprob, ps)
  }
  nodeprob <- normalizeProbSlot(node, parents, categories,
                                seq(1,length(parents)), nodeprob)
  return(nodeprob)
}

setMethod("cnNodeLoglik", c("catNetwork"), 
          function(object, node, data, perturbations=NULL) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of categories")

            if(is.data.frame(data)) {
              data <- as.matrix(t(data))
              if(!is.null(perturbations)) { 
                if(!is.data.frame(perturbations))
                  stop("Perturbations should be a data frame")
                perturbations <- as.matrix(t(perturbations))
              }
            }

            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
              stop("The number of nodes in  the object and data should be equal")
            
            rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows/cols should be named after the nodes of the object.")

            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              ## keep the matrix format !!
              data <- as.matrix(data[norder,])
              if(!is.null(perturbations))
                perturbations <- as.matrix(perturbations[norder,])
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
              if(is.numeric(node))
                node <- sapply(node, function(nn) return(which(norder == nn)))
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The data names should correspond to the object nodes.")

            if(is.character(node))
              node <- sapply(node, function(nn) return(which(rownames == nn)))
            if(!is.numeric(node) || sum(node < 1) || sum(node > object@numnodes))
              stop("Wrong node")
            node <- as.integer(node)
 
            r <- .categorizeSample(data, perturbations, object, ask=FALSE)
            data <- r$data
            perturbations <- r$perturbations
            if(object@maxCategories < r$maxCategories)
              stop("Data has more categories than the object")
             
            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]
            
            if(TRUE) {
              loglik <- .Call("ccnNodeLoglik", 
                              object, node, data, perturbations, 
                              PACKAGE="catnet")
              return(loglik)
            }
            else
              return(.nodeLikelihood(node,
                                     object@parents[[node]], object@categories, 
                                     seq(1,length(object@parents[[node]])),
                                     object@probabilities[[node]], data))
          })

setMethod("cnNodeLoglikError", c("catNetwork"), 
          function(object, node, data, perturbations=NULL) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of categories")

            if(is.data.frame(data))
              data <- as.matrix(t(data))

            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
              stop("The number of nodes in  the object and data should be equal")
            
            rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows/cols should be named after the nodes of the object.")
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              ## keep the matrix format !!
              data <- as.matrix(data[norder,])
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The data names should correspond to the object nodes.")

            if(!is.numeric(node) || node < 1 || node > object@nodes)
              stop("Wrong node")
            node <- as.integer(node)
              
            r <- .categorizeSample(data, NULL, object, ask=FALSE)
            data <- r$data
            if(object@maxCategories < r$maxCategories)
              stop("Data has more categories than the object")
            
            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]
            
            loglik <- .Call("ccnNodeLoglikError", 
                            object, node, data, NULL, 
                            PACKAGE="catnet")
            return(loglik)
          })
