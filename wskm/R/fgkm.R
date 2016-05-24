fgkm <- function(x, centers, groups, lambda, eta, maxiter=100, delta=0.000001, maxrestart=10,seed=-1) 
{
  if (missing(centers))
    stop("the number or initial clusters 'centers' must be provided")

  if(seed<=0){
    seed <-runif(1,0,10000000)[1]
  }

  vars <- colnames(x)
  
  nr <-nrow(x) # nrow() return a integer type
  nc <-ncol(x) # integer

  if (is.data.frame(centers) || is.matrix(centers))
  {
    init <- TRUE
    k <- nrow(centers)
  }
  else
  {
    init <- FALSE
    k <- centers
    centers <- double(k * nc)
  }
  
  # get the setting of feature group
  if (is.character(groups) && length(groups) == 1) {
    G <- .C("parseGroup",as.character(groups),numGroups=integer(1), groupInfo=integer(nc),PACKAGE="wskm")
  } else if (is.vector(groups) && length(groups) == nc) {
    G <- list()
    grps <- as.factor(groups)
    groupNames <- levels(grps)
    G$numGroups <- nlevels(grps)
    G$groupInfo <- as.integer(as.integer(grps) - 1)
  }

  set.seed(seed)
  Z <- .C("fgkm",
          x = as.double(as.matrix(x)),
          nr,
          nc,
          k = as.integer(k),
          lambda = as.double(lambda),
          eta = as.double(eta),
          G$numGroups,
          G$groupInfo,
          delta = as.double(delta),
          maxIterations = as.integer(maxiter),
          maxRestarts = as.integer(maxrestart),
          as.logical(init),
#          seed,
          cluster = integer(nr),
          centers=as.double(as.matrix(centers)),
          featureWeight = double(k * nc),
          groupWeight = double(k * G$numGroups),
          iterations = integer(1),
          restarts = integer(1),
          totiters = integer(1),
          totalCost = double(1),
          totss = double(1),
		  withiness = double(k),
          PACKAGE="wskm"
          )
       
  centers <- matrix( Z$centers)
  dim(centers) <- c(k, nc)
  colnames(centers) <- vars
  
  featureWeight <- matrix(Z$featureWeight)
  dim(featureWeight) <- c(k, nc)
  colnames(featureWeight) <- vars
  
  groupWeight <- matrix(Z$groupWeight)
  dim(groupWeight) <- c(k, G$numGroups )
  colnames(groupWeight) <- 1:ncol(groupWeight)
  
  ignore <- which(rowSums(centers==0) == ncol(centers))
  if (length(ignore)) {
    centers <- centers[-ignore,, drop=FALSE]
    featureWeight <- featureWeight[-ignore,, drop=FALSE]
  }
  
  rownames(centers) <- 1:nrow(centers)
  rownames(featureWeight) <- 1:nrow(featureWeight)
  rownames(groupWeight) <- 1:nrow(groupWeight)
  
  cluster <- Z$cluster + 1
  
  size <- aggregate(cluster, list(cluster=cluster), length)[[2]]
  
  result <- list(cluster = cluster,
                 centers = Z$centers,
                 totss = Z$totss, 
                 withinss = Z$withinss, 
                 tot.withinss = sum(Z$withiness), 
                 betweenss = Z$totss-sum(Z$withinss),
                 size = size,
                 iterations = Z$iterations,
                 restarts = Z$restarts,
                 totiters=Z$totiters,
                 featureWeight = Z$featureWeight,
                 groupWeight = Z$groupWeight)
  
  dim(result$centers) <- c(k, nc)
  dim(result$featureWeight) <- c(k, nc)
  dim(result$groupWeight) <- c(k, G$numGroups)
  
  class(result) <- c("kmeans", "fgkm")
  return(result)
}
