########################################################################
# Categorical Network Class Methods
# Searching

cnGenOrder <- function(order, shuffles, bjump = FALSE) {
  if(shuffles > 0) {
    neworder <- order
    nn <- length(order)
    for(sh in 1:shuffles) {
      order <- neworder
      n1 <- floor(1+runif(1,0,1)*nn)
      if(bjump) {
      	n2 <- n1
        while(n2==n1)
          n2 <- floor(1+runif(1,0,1)*nn)
      }
      else {
        n2 <- n1 + 1
        if(n1 == nn)
          n2 <- 1
      }
      if(n1 < n2) {
        if(n1 > 1)
          neworder[1:(n1-1)] <- order[1:(n1-1)]
        neworder[n1:(n2-1)] <- order[(n1+1):n2]
        neworder[n2] <- order[n1]
        if(n2 < nn)
          neworder[(n2+1):nn] <- order[(n2+1):nn]
      }
      else {
        if(n2 > 1)
          neworder[1:(n2-1)] <- order[1:(n2-1)]
        neworder[(n2+1):n1] <- order[n2:(n1-1)]
        neworder[n2] <- order[n1]
        if(n1 < nn)
          neworder[(n1+1):nn] <- order[(n1+1):nn]
      }
    }
  }
  else{ 
    neworder <- sample(1:length(order))
  }
  neworder
}

cnSearchOrder <- function(data, perturbations = NULL,
                          maxParentSet = 0, parentSizes = NULL, 
                          maxComplexity = 0,
                          nodeOrder = NULL, nodeCats = NULL, 
                          parentsPool = NULL, fixedParents = NULL,
                          edgeProb = NULL, 
                          echo = FALSE) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  fast <- TRUE
  
  t1 <- proc.time()
  
  if(is.matrix(data)) {
    numnodes <- nrow(data)
    numsamples <- ncol(data)
  }
  else {
    numnodes <- ncol(data)
    numsamples <- nrow(data)
  }
  if(numnodes < 1 || numsamples < 1)
    stop("Invalid sample")

  if(!is.null(nodeCats))
    if(!is.list(nodeCats) || length(nodeCats) != numnodes)
      stop("Wrong nodeCats parameter")

  r <- .categorizeSample(data, perturbations, object=NULL, nodeCats=nodeCats, ask=TRUE)
  data <- r$data
  perturbations <- r$perturbations
  categories <- r$categories
  maxCategories <- r$maxCategories

  nodenames <- rownames(data)    
  if(length(nodenames) < numnodes)
    nodenames <- seq(1,numnodes)
  if(length(nodenames) > 1)
    for(i in 2:length(nodenames))
      if(length(which(nodenames[1:(i-1)] == nodenames[i])) > 0) {
        stop("Repeated node names ", nodenames[i])
      }
  
  maxParentSet <- as.integer(maxParentSet)
  if(maxParentSet < 1) {
    if(!is.null(parentSizes))
      maxParentSet <- as.integer(max(parentSizes))
    if(maxParentSet < 1) 
      maxParentSet <- 1
  }
  
  if(!is.null(parentSizes) && length(parentSizes) == numnodes) {
    parentSizes <- as.integer(parentSizes)
    parentSizes[parentSizes<0] <- 0
    parentSizes[parentSizes>maxParentSet] <- maxParentSet
  }
  else
    parentSizes <- NULL

  if(!is.null(edgeProb)) {
    if(!is.matrix(edgeProb) || nrow(edgeProb) != numnodes || ncol(edgeProb) != numnodes)
      stop("edgeProb should be square matrix of length the number of nodes")
    for(i in 1:numnodes) {
      edgeProb[i, edgeProb[i,] < 0] <- 0
      edgeProb[i, edgeProb[i,] > 1] <- 1
    }

    bones <- 0
    for(i in 1:numnodes) {
      for(j in 1:numnodes) {
        if(is.na(edgeProb[i,j]) || edgeProb[i,j] < 0 || edgeProb[i,j] > 1)
           edgeProb[i,j] <- 0.5
        if(edgeProb[i,j] == 1)
          bones <- bones + 1
      }
    }
    if(bones) {
      if(is.null(fixedParents))
        fixedParents <- vector("list", numnodes)
      for(i in 1:numnodes) {
        for(j in 1:numnodes) {
          if(edgeProb[i,j] == 1 && length(fixedParents[[i]]) < maxParentSet) {
            fixedParents[[i]] <- c(fixedParents[[i]], j)
            ##cat("fixed ", i, ": ", fixedParents[[i]], "\n")
          }
        }
      }
    }
  }

  if(!is.null(parentsPool)) {
    if(!is.list(parentsPool) || length(parentsPool) != numnodes)
      stop("Wrong parentsPool parameter")
    for(i in 1:numnodes)
      if(is.character(parentsPool[[i]])) {
        ids <- NULL
	for(ch in parentsPool[[i]]) {
          id <- which(nodenames == ch)
          if(length(id) <= 0)
            stop("Invalid parentsPool")
          ids <- c(ids, id[1])
        }
        parentsPool[[i]] <- as.integer(ids)
      }
  }
    
  if(!is.null(fixedParents)) {
    if(!is.list(fixedParents) || length(fixedParents) != numnodes)
      stop("Wrong fixedParents parameter")
    for(i in 1:numnodes) {
      if(length(fixedParents[[i]]) > maxParentSet) {
        fixedParents[[i]] <- fixedParents[[i]][1:maxParentSet]
        warning("fixedParents for node ", i, " exceeds the maximum parent size. It's reduced to ", maxParentSet)
      }
      if(is.character(fixedParents[[i]])) {
        ids <- NULL
	for(ch in fixedParents[[i]]) {
          id <- which(nodenames == ch)
          if(length(id) <= 0)
            stop("Invalid fixedParents")
          ids <- c(ids, id[1])
        }
        fixedParents[[i]] <- as.integer(ids)
      }
    }
  }

  if(is.null(nodeOrder))
    nodeOrder <- 1:numnodes
  if(!is.null(nodeOrder)) {
    if(length(nodeOrder) != numnodes)
      stop("Invalid nodeOrder")
    if(is.character(nodeOrder)) {
      nodeOrder <- sapply(nodeOrder, function(str) {
        id <- which(nodenames == str)
        if(length(id) <= 0)
          stop("Invalid nodeOrder")
        return(id[1])
        })
    }
  }

  if(maxComplexity <= 0)
    maxComplexity <- as.integer(numnodes * exp(log(maxCategories)*maxParentSet) * (maxCategories-1))
  minComplexity <- sum(sapply(categories, function(cat) (length(cat)-1)))
  if(maxComplexity < minComplexity) {
    cat("set maxComplexity to ", minComplexity, "\n")
    maxComplexity <- minComplexity
  }
  maxComplexity <- as.integer(maxComplexity)
  
  catIndices <- NULL
  if(!is.null(nodeCats)) {
    catIndices <- lapply(1:numnodes, function(i) 1:length(categories[[i]]))
  }
    
  .Call("ccnReleaseCache", PACKAGE="catnet")
  bestnets <- .Call("ccnOptimalNetsForOrder", 
                    data, perturbations, 
                    as.integer(maxParentSet), as.integer(parentSizes),
                    as.integer(maxComplexity),
                    as.integer(nodeOrder), catIndices, 
                    parentsPool, fixedParents, edgeProb, 
                    ## no cache
                    FALSE, 
                    echo, 
                    PACKAGE="catnet")
  if(length(nodenames) == numnodes && length(bestnets) > 0) {
    for(i in 1:length(bestnets)) {
      if(is.null(bestnets[[i]])) {
        warning("No network")
        next
      }
      bestnets[[i]]@nodes <- nodenames[nodeOrder]
    }
  }
  
  if(echo)
    cat("Collating ", length(bestnets), " networks...\n")
  
  eval <- new("catNetworkEvaluate", numnodes, numsamples, length(bestnets))
  
  for(i in 1:length(bestnets)) {
    if(is.null(bestnets[[i]]))
      next
    eval@complexity[i] <- bestnets[[i]]@complexity
    eval@loglik[i] <- bestnets[[i]]@likelihood
    ## bestnets are ordered according to `nodidx'
    ## must set their categorical values
    bestnets[[i]]@categories <- categories[nodeOrder]

    ## reorder eval@nets[[nn]]'s nodes to match data's nodes
    enetnodes <- bestnets[[i]]@nodes
    ##cat("nodenames = ", nodenames,"\n")
    ##cat("enetnodes = ", enetnodes,"\n")
    if(length(nodenames) == numnodes) {
      ord <- sapply(nodenames, function(c) {
        id <- which(enetnodes==c)
        if(length(id)>0)
          return(id[1])
	stop("nodes do not match")
        })
      if(sum(ord != 1:numnodes) > 0) {
        bestnets[[i]] <- cnReorderNodes(bestnets[[i]], ord)
      }
    }
  }
  eval@nets <- bestnets

  t2 <- proc.time()
  eval@time <- eval@time + as.numeric(t2[3] - t1[3])
  
  return(eval)
}

.searchSA <- function(eval, 
                      numnodes, numsamples, 
                      data, perturbations, 
                      categories, maxCategories,
                      maxParentSet, parentSizes, 
                      maxComplexity,
                      startorder,
                      catIndices, 
                      parentsPool, fixedParents, maxParentsPool, 
                      edgeProb, dirProb, 
                      selectMode, 
                      tempStart, tempCoolFact, tempCheckOrders, 
                      maxIter, orderShuffles, stopDiff, 
                      numThreads, echo) {
  
  nodenames <- rownames(data)
  .Call("ccnReleaseCache", PACKAGE="catnet")
  optnets <- .Call("ccnOptimalNetsSA",
                   nodenames, 
                   data, perturbations, 
                   as.integer(maxParentSet), as.integer(parentSizes),
                   as.integer(maxComplexity),
                   catIndices, 
                   parentsPool, fixedParents, as.integer(maxParentsPool), 
                   edgeProb, dirProb, 
                   selectMode, as.integer(startorder),
                   tempStart, tempCoolFact, as.integer(tempCheckOrders), 
                   as.integer(maxIter), orderShuffles, stopDiff, 
                   ## threads
                   as.integer(numThreads),
                   ## cache
                   TRUE, 
                   echo, 
                   PACKAGE="catnet")
  
  if(echo)
    cat("Collating ", length(optnets), " networks...\n")
  
  for(nn in 1:length(optnets)) {
    if(is.null(optnets[[nn]])) {
      warning("No network")
      next
    }
    enetnodes <- optnets[[nn]]@nodes
    ord <- sapply(nodenames, function(c) which(enetnodes==c))
    if(length(ord) != numnodes)
      stop("wrong node names: ", enetnodes)
    if(sum(ord != 1:numnodes) > 0)    
      optnets[[nn]] <- cnReorderNodes(optnets[[nn]], ord)
    ## now the nodes in the network are ordered as in the data
    optnets[[nn]]@categories <- categories
  }
  
  ## replace existing best networks    
  if(!is.null(eval) && length(eval@nets) > 0) {
    enets <- optnets
    for(nn in 1:length(enets)) {
      if(is.null(enets[[nn]]) || !is(enets[[nn]], "catNetwork"))
        next
      bnet <- cnFind(eval@nets, enets[[nn]]@complexity)
      if(!is.null(bnet) && bnet@complexity == enets[[nn]]@complexity) 
        if(bnet@likelihood > enets[[nn]]@likelihood)
          ##replace with better
          enets[[nn]] <- bnet
    }
    eval@nets <- enets
  }
  else {
    eval@nets <- optnets
  }

  eval@complexity <- sapply(eval@nets, function(net) net@complexity)
  eval@loglik <- sapply(eval@nets, function(net) net@likelihood)
  
  return(eval)
}

cnSearchSA <- function(data, perturbations=NULL,
                       maxParentSet=0, parentSizes = NULL, 
                       maxComplexity = 0, nodeCats = NULL, 
                       parentsPool = NULL, fixedParents = NULL, 
                       edgeProb = NULL, dirProb = NULL, 
                       selectMode = "BIC", 
                       tempStart = 1, tempCoolFact = 0.9, tempCheckOrders = 10, 
                       maxIter = 100, orderShuffles = 1, stopDiff = 0, 
                       numThreads = 2, 
                       priorSearch = NULL,  ## catNetworkEvaluate
                       echo=FALSE) {

  t1 <- proc.time()
  
  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  if(is.matrix(data)) {
    numnodes <- nrow(data)
    numsamples <- ncol(data)
  }
  else {
    numnodes <- ncol(data)
    numsamples <- nrow(data)
  }

  if(numnodes < 1 || numsamples < 1)
    stop("Insufficient data")

  if(!is.null(nodeCats))
    if(!is.list(nodeCats) || length(nodeCats) != numnodes)
      stop("Wrong nodeCats parameter")

  r <- .categorizeSample(data, perturbations, object=NULL, nodeCats=nodeCats, ask=TRUE)
  data <- r$data
  perturbations <- r$perturbations
  categories <- r$categories
  maxCategories <- r$maxCategories

  nodenames <- rownames(data)
  if(length(nodenames) < numnodes)
    nodenames <- seq(1, numnodes)
  if(length(nodenames) > 1)
    for(i in 2:length(nodenames))
      if(length(which(nodenames[1:(i-1)] == nodenames[i])) > 0) {
        stop("Repeated node names ", nodenames[i])
      }

  maxParentSet <- as.integer(maxParentSet)
  if(maxParentSet < 1) {
    if(!is.null(parentSizes))
      maxParentSet <- as.integer(max(parentSizes))
    if(maxParentSet < 1) 
      maxParentSet <- 1
  }

  if(!is.null(parentSizes)) {
    parentSizes <- as.integer(parentSizes)
    parentSizes[parentSizes<0] <- 0
    parentSizes[parentSizes>maxParentSet] <- maxParentSet
  }
  else
    parentSizes <- NULL

  if(!is.null(parentsPool)) {
    if(!is.list(parentsPool) || length(parentsPool) != numnodes)
      stop("Wrong parentsPool parameter")
    for(i in 1:numnodes)
      if(is.character(parentsPool[[i]])) {
        ids <- NULL
	for(ch in parentsPool[[i]]) {
          id <- which(nodenames == ch)
          if(length(id) <= 0)
            stop("Invalid parentsPool")
          ids <- c(ids, id[1])
        }
        parentsPool[[i]] <- as.integer(ids)
      }
  }

  if(!is.null(fixedParents)) {
    if(!is.list(fixedParents) || length(fixedParents) != numnodes)
      stop("Wrong fixedParents parameter")
    for(i in 1:numnodes) {
      if(length(fixedParents[[i]]) > maxParentSet) {
        fixedParents[[i]] <- fixedParents[[i]][1:maxParentSet]
        warning("fixedParents for node ", i, " exceeds the maximum parent size. It's reduced to ", maxParentSet)
      }
      if(is.character(fixedParents[[i]])) {
        ids <- NULL
	for(ch in fixedParents[[i]]) {
          id <- which(nodenames == ch)
          if(length(id) <= 0)
            stop("Invalid fixedParents")
          ids <- c(ids, id[1])
        }
        fixedParents[[i]] <- as.integer(ids)
      }
    }
  }
  
  if(!is.null(priorSearch)) {
    if(!is(priorSearch, "catNetworkEvaluate")) 
      stop("'priorSearch' should be a valid catNetworkEvaluate object or NULL")
    if(priorSearch@numnodes != numnodes || priorSearch@numsamples != numsamples)
      stop("priorSearch's number of nodes and sample size are not compatible with those of the data")
  }
  
  if(tempStart < 0) {
    warning("tempStart is set to 0")
    tempStart <- 0
  }
  if(tempCoolFact < 0) {
    warning("tempCoolFact is set to 0")
    tempCoolFact <- 0
  }
  if(tempCoolFact > 1) {
    warning("tempCoolFact is set to 1")
    tempCoolFact <- 1
  }
  if(tempCheckOrders < 1) {
    warning("tempCheckOrders is set to 1")
    tempCheckOrders <- 1
  }
  if(maxIter < tempCheckOrders) {
    warning("maxIter is set to tempCheckOrders")
    maxIter <- tempCheckOrders
  }

  catIndices <- NULL
  if(!is.null(nodeCats)) {
    catIndices <- lapply(1:numnodes, function(i) 1:length(categories[[i]]))
  }
  
  if(FALSE && is.null(dirProb) && !is.null(perturbations)) {
    mat <- .Call("ccnKLPairwise", 
                  data, perturbations, 
                  PACKAGE="catnet")
    klmat <- matrix(mat, nrow=numnodes, ncol=numnodes)
    dirProb <-t(klmat)/(klmat+t(klmat))
    rownames(dirProb)<-nodenames
    colnames(dirProb)<-nodenames
  }

  if(!is.null(edgeProb)) {
    if(length(colnames(edgeProb)) != numnodes || length(rownames(edgeProb)) != numnodes)
      stop("edgeProb matrix should be with named columns and rows")
    edgeProb <- edgeProb[nodenames, nodenames]
  }

  if(!is.null(dirProb)) {
    if(length(colnames(dirProb)) != numnodes || length(rownames(dirProb)) != numnodes)
      stop("dirProb matrix should be with named columns and rows")
    dirProb <- dirProb[nodenames, nodenames]
    dirProb <- dirProb / (dirProb+t(dirProb))
    for(i in 1:numnodes) {
      for(j in 1:numnodes) {
        if(is.na(dirProb[i,j]) || is.na(dirProb[i,j]) || dirProb[i,j] < 0 || dirProb[i,j] > 1)
           dirProb[i,j] <- 0.5
        if(dirProb[i,j] > 1-1e-8)
           dirProb[i,j] <- 1-1e-8
        if(dirProb[i,j] < 1e-8)
           dirProb[i,j] <- 1e-8
      }
    }
  }
  
  if(!is.null(fixedParents)) {
    for(i in 1:numnodes)
      if(length(fixedParents[[i]]) > maxParentSet) {
        fixedParents[[i]] <- fixedParents[[i]][1:maxParentSet]
        warning("fixedParents for node ", i, " exceeds the maximum parent size. It's reduced to ", maxParentSet)
      }
  }
  
  if(maxComplexity <= 0)
    maxComplexity <- as.integer(numnodes * exp(log(maxCategories)*maxParentSet) * (maxCategories-1))
  minComplexity <- sum(sapply(categories, function(cat) (length(cat)-1)))
  if(maxComplexity < minComplexity) {
    warning("set maxComplexity to ", minComplexity)
    maxComplexity <- minComplexity
  }

  maxParentsPool <- 0
  if(maxParentsPool > 0) {
    maxParentsPool <- as.integer(maxParentsPool)
    if(maxParentsPool < maxParentSet)
      maxParentsPool <- maxParentSet
    if(maxParentsPool > numnodes)
      maxParentsPool <- numnodes
  }
  
  if(is.null(priorSearch) || length(priorSearch@nets) < 1) {
    optorder <- sample(1:numnodes)
    eval <- new("catNetworkEvaluate", numnodes, numsamples, 0)
    eval@time <- 0
  }
  else {
    eval <- priorSearch
    optnets <- priorSearch@nets
    if(length(optnets) < 1)
      stop("'priorSearch' has no networks")

    if(selectMode == "AIC") {
      optnet <- cnFindAIC(optnets, numsamples)
    }
    else {
      if(selectMode == "BIC") {
        optnet <- cnFindBIC(optnets, numsamples)
      }
      else {
        optnet <- cnFind(optnets, maxComplexity)
      }
    }
    if(is.null(optnet))
      optnet <- optnets[[1]]
    if(is.null(optnet)) {
      optorder <- sample(1:numnodes)
      warning("no valid network found - priorSearch ignored")
    }
    optorder <- cnOrder(optnet)

    ## make sure the optnet comes is data compatible
    if(prod(nodenames == optnet@nodes) == 0)
      stop("The data names should correspond to the prior search nodes.")
  }
  
  eval <- .searchSA(eval,
                    as.integer(numnodes), as.integer(numsamples),
                    data, perturbations, 
                    categories, as.integer(maxCategories),
                    as.integer(maxParentSet), parentSizes, 
                    as.integer(maxComplexity),
                    as.integer(optorder),
                    catIndices, 
                    parentsPool, fixedParents, as.integer(maxParentsPool), 
                    edgeProb, dirProb, 
                    selectMode, 
                    tempStart, tempCoolFact, as.integer(tempCheckOrders), 
                    as.integer(maxIter), orderShuffles, stopDiff, 
                    as.integer(numThreads), echo)
  
  t2 <- proc.time()
  eval@time <- eval@time + as.numeric(t2[3] - t1[3])

  return(eval)
}
