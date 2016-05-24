########################################################################
# Categorical Network Class Methods
# Histograms

##setMethod("cnParHist", "list", 
##          function(objectlist) {
##            if(!is(objectlist, "list"))
##              stop("A list of catNetworks should be specified.")
##            if(length(objectlist)==0 || !is(objectlist[[1]], "catNetwork"))
##              stop("At least one catNetworks should be specified.")
##            return(parHisto(objectlist))
##          })

parHisto <- function(objectlist, norder = NULL) {
  if(is(objectlist, "catNetwork")) {
    n <- objectlist@numnodes
  if(is.null(norder))
    norder <- seq(1, n)
    return(matParents(objectlist, norder))
  }
  
  n <- objectlist[[1]]@numnodes
  if(is.null(norder))
    norder <- seq(1, n)
  
  i <- 1
  nnets <- length(objectlist)
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    if(object@numnodes != n)
      stop("Networks should have the same number of nodes.")
    mpar <- matParents(object, norder)
    if(i==1)
      mhisto <- mpar
    else
      mhisto <- mhisto + mpar

    i <- i + 1
  }
 
  return(mhisto)
}
 
cnSearchHist <- function(data, pert=NULL,  
                         maxParentSet=1, parentSizes = NULL,
                         maxComplexity=0, nodeCats = NULL,  
                         parentsPool = NULL, fixedParents = NULL,
                         score = "BIC", weight="loglik", 
                         maxIter = 32, numThreads = 2, echo=FALSE) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")
  
  if(is.matrix(data)) {
    numnodes <- nrow(data)
    numsamples <- ncol(data)
    nodenames <- rownames(data)
  }
  else {
    numnodes <- ncol(data)
    numsamples <- nrow(data)
    nodenames <- colnames(data)
  }
  
  if(numnodes < 1 || numsamples < 1)
    stop("No valid sample is specified.")

  if(length(nodenames) < numnodes) {
    nodenames <- seq(1, numnodes)
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
  
  r <- .categorizeSample(data, pert, object=NULL, nodeCats=nodeCats, ask=TRUE)
  data <- r$data
  pert <- r$pert
  cats <- r$cats
  maxcats <- r$maxcats

  catIndices <- NULL
  if(!is.null(nodeCats)) {
    catIndices <- lapply(1:numnodes, function(i) 1:length(cats[[i]]))
  }
  
  if(maxComplexity <= 0)
    maxComplexity <- as.integer(numnodes * exp(log(maxcats)*maxParentSet) * (maxcats-1))
  minComplexity <- sum(sapply(cats, function(cat) (length(cat)-1)))
  if(maxComplexity < minComplexity)
    maxComplexity <- minComplexity
  
  numThreads <- as.integer(numThreads)
  if(numThreads < 1)
    numThreads <- 1

  maxIter <- as.integer(maxIter)
  if(maxIter < numThreads)
    maxIter <- numThreads

  nweight <- 0
  if(weight=="loglik")
    nweight <- 1
  if(weight=="score")
    nweight <- 2
  
  ## call the C-function
  .Call("ccnReleaseCache", PACKAGE="sdnet")
  vhisto <- .Call("ccnParHistogram", 
                  data, pert, 
                  as.integer(maxParentSet), as.integer(parentSizes),
                  as.integer(maxComplexity),
                  catIndices, 
                  parentsPool, fixedParents,
                  score, nweight, as.integer(maxIter),
                  as.integer(numThreads), 
                  ## cache
                  TRUE, 
                  echo, 
                  PACKAGE="sdnet")
  
  mhisto <- matrix(vhisto, numnodes, numnodes)
  rownames(mhisto)<-nodenames
  colnames(mhisto)<-nodenames

  return(mhisto)
}
