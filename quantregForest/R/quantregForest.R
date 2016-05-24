"quantregForest" <-
function(x,y, ...){

  ## Some checks 
  if(! class(y) %in% c("numeric","integer") )
    stop(" y must be numeric ")
  
  if(is.null(nrow(x)) || is.null(ncol(x)))
    stop(" x contains no data ")
    

  if( nrow(x) != length(y) )
    stop(" predictor variables and response variable must contain the same number of samples ")

  if (any(is.na(x))) stop("NA not permitted in predictors")
  if (any(is.na(y))) stop("NA not permitted in response")

  
  
  ## Check for categorial predictors with too many categories (copied from randomForest package)
   if (is.data.frame(x)) {
        ncat <- sapply(x, function(x) if(is.factor(x) && !is.ordered(x))
                       length(levels(x)) else 1)
      } else {
        ncat <- 1
    }
    maxcat <- max(ncat)
    if (maxcat > 32)
        stop("Can not handle categorical predictors with more than 32 categories.")

  
  ## Note that crucial parts of the computation
  ## are only invoked by the predict method
  cl <- match.call()
  cl[[1]] <- as.name("quantregForest")

  qrf <- randomForest( x=x,y=y ,...)
  nodesX <- attr(predict(qrf,x,nodes=TRUE),"nodes")
  rownames(nodesX) <- NULL
  nnodes <- max(nodesX)
  ntree <- ncol(nodesX)
  n <- nrow(x)
  valuesNodes <- matrix(nrow=nnodes,ncol=ntree)
  for (tree in 1:ntree){
      shuffledNodes <- nodesX[rank(ind <- sample(1:n,n)),tree]
      useNodes <- sort(unique(as.numeric(shuffledNodes)))
      valuesNodes[useNodes,tree] <- y[ind[match(useNodes,shuffledNodes )]]
  }
  
  class(qrf) <- c("quantregForest","randomForest")

  qrf[["call"]] <- cl
  qrf[["valuesNodes"]] <- valuesNodes

  
  return(qrf)
}
