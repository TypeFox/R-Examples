### subfunction for final computation of within cluster projected energies and sparsity coefficient
# compute projected new subspace and projected energy merge of two clusters
compute.final.projen <- function(i, x, act.clustering, allpts = FALSE){
  if (!allpts) { 
    clusteri    <- x[which(act.clustering$cluster == i),] # cluster i
    centpro     <- t(apply(clusteri %*% act.clustering$subspaces[[i]], 1, function(z, zen) return(z-zen), zen = t(act.clustering$centers[i,])%*%act.clustering$subspaces[[i]])) # cluster in subspace shifted by its mean 
    projen      <- sum(apply(centpro, 1, function(z) return(sum(z^2)))) / nrow(centpro)  # projected energy: average distances from center in subspace
    }
  if (allpts) { 
    cntr        <- apply(x, 2, mean)
    centpro     <- t(apply(x %*% act.clustering$subspaces[[i]], 1, function(z, zen) return(z-zen), zen = t(cntr)%*%act.clustering$subspaces[[i]])) # cluster in subspace shifted by its mean 
    projen      <- sum(apply(centpro, 1, function(z) return(sum(z^2)))) / nrow(centpro)  # projected energy: average distances from center in subspace
    }    
  return(projen)
  }


################################
### definition of orclus method:
################################

orclus <- function (x, ...) 
    UseMethod("orclus")


###################################
### main orclus.default - function:
###################################

orclus.default <- function(x, k, l, k0, a = 0.5, inner.loops = 1, verbose = TRUE, ...){
  # x:            data input matrix
  # k:            number of clusters
  # l:            final subspace dimension
  # k0:           number of initial clusters
  # a:            reduction rate for subspace dimension over the iterations
  # inner.loops:  number of iterations for each pair of 'number of clusters' and 'subspace diension' during the iteration process
  # verbose:      logical specifying whether iteration process should be printed
  
  if (is.data.frame(x)) x <- as.matrix(x)
  
  # original dimensions
  d <- ncol(x)
  n <- nrow(x)
  
  # some preliminary error checks
  if (d < l) stop("Subspace dimension should be specified lower than the input dimension!")
  if (d == l){warning("Specified subspace dimension equal to original original dimension. kmeans() applied instead.")
              return(kmeans(x, centers = k))}
  if (d == 1) stop("For dimension of one no subspaces can be computed!")
  if (l < 1) stop("Subspace dimension l should be at least 1!")
  if (n == l) stop("No clustering is necessary if data consists of only one observation!")
  if (n < k*l) stop("Data must consist of more than k*l observations to compute dimensionality reduction for the clusters!")
  if (n < k) stop("Ther can't be more clusters than objects!")
  if (k == 1) stop("Clustering is not meaningful for only one cluster! For dimensionality reduction use principal component analysis.")
  if (!is.numeric(x)) stop("All input variables must be numeric!")
  if (a > 1 | a < 0) stop("Parameter a must be from the interval (0,1)!")
  if (k0 <= k) stop("For the initial number of clusters k0 > k is supposed!")
  
  # reduction factor for the number of clusters
  b <- exp(-log(d/l)*log(1/a)/log(k0/k))
  niters <- 1 + ceiling(log(k/k0, a)) 
  
  # initialization 
  act.d <- d  # current dimension and number of clusters
  act.k <- k0
  kc    <- k0 # initialization of new number of clusters and subspace dimension for the first iteration to k and l 
  lc    <- d

  act.clustering <- list()
  
  # begin iteration 
  for(iter in 1:niters){ 
    if (iter == 1 | act.d > l | act.k > k) {
       for(ilp in 1:inner.loops){ 
         if (verbose & ilp == 1) cat("iteration         :",iter,"\n")
              
         #1) assign clusters:
         if (iter == 1 & ilp == 1) {     
            if (verbose) cat("Initialization with", act.k, "clusters \n")
            clusinit                <- kmeans(x, max(k,round(kc)))
            act.clustering$cluster  <- clusinit$cluster
            act.clustering$centers  <- clusinit$centers
            act.clustering$size     <- clusinit$size 
            remove(clusinit)
            }
         if (iter > 1){
            clusupdate <- reassign(x, act.clustering)
            act.clustering$cluster   <- clusupdate$cluster
            act.clustering$centers   <- clusupdate$centers
            act.clustering$size      <- clusupdate$size
            act.clustering$subspaces <- clusupdate$subspaces
            remove(clusupdate)         
            }
             
         #if (verbose & ilp == 1) cat("current ratio (dimension) :", lc, "\n") 
         
         #2) find new subspaces
         act.clustering$subspaces <- findvectors(x, act.clustering, dimen = max(l,round(lc)))
  
         if(is.vector(act.clustering$subspaces[[1]])) act.d <- 1
         if(is.matrix(act.clustering$subspaces[[1]])) act.d <- dim(act.clustering$subspaces[[1]])[2]
         if (verbose & ilp == 1) cat("Actual Subspace dimension :", act.d, "\n") 
         #if (verbose) cat("Loop                      :", ilp,"\n")
         } 
        
       #3) set new dimension and cluster number for next iteration
       knew <- kc * a
       lnew <- lc * b
       
       #4) merge clusters
       #if (verbose) cat("current ratio (cluster)   :", knew, "\n") 
       act.clustering <- clmerge(x, act.clustering, max(k, round(knew)))     

       act.k <- length(table(act.clustering$cluster))
       kc    <- knew
       lc    <- lnew
       if (verbose) cat("New number of clusters    :", act.k, "\n\n")
       }
  }
  # final reassignment of the clusters
  if (verbose) cat("Final reassigment...\n")   
  clusupdate <- reassign(x, act.clustering)
  act.clustering$cluster   <- clusupdate$cluster
  act.clustering$centers   <- clusupdate$centers
  act.clustering$size      <- clusupdate$size
  act.clustering$subspaces <- clusupdate$subspaces
  remove(clusupdate) 
  
  # check if due to removal of empty clusters the final cluster size is lower than the prespecified k
  if (length(act.clustering$size) < k){ 
    k <- length(act.clustering$size)
    warning("Due to empty classed after merge less than k classes result!")
    cat("New cluster number due to empty classes:", k, "\n\n")
    }
        
  # compute projected energy per cluster and ad it to result list ...this corresponds to withinss from kmeans
  within.projens <- sapply(1:k, compute.final.projen, x = x, act.clustering = act.clustering, allpts = FALSE)

  # compute sparsity coefficient and add it to result list 
  projens.wholedata     <- sapply(1:k, compute.final.projen, x = x, act.clustering = act.clustering, allpts = TRUE)
  sparsity.coefficient  <- mean(within.projens / projens.wholedata)

  act.clustering$subspace.dimension   <- l
  act.clustering$within.projens       <- within.projens
  act.clustering$sparsity.coefficient <- sparsity.coefficient
  act.clustering$orclus.call <- match.call()
  act.clustering$orclus.call[[1]] <- as.name("orclus")
  
  # create result-object of class orclus
  class(act.clustering) <- "orclus"
  return(act.clustering)   
  }



########################
### predict function
########################

predict.orclus <- function(object, newdata, ...){
  # several error checks 
  if (class(object) != "orclus") stop("Object must be of class oclus!")
  if(is.data.frame(newdata)) newdata <- as.matrix(newdata)
  if(!is.matrix(newdata)) stop("Newdata must be of type matrix or data frame!")
  if(dim(object$centers)[2] != dim(newdata)[2]) stop("Newdata must be of same dimension as the explanatory input variables data set!")
  
  clids         <- 1:length(table(object$cluster))

  # compute distances of the observations to each clusterin the correspondin subspaces
  comp.dists         <- function(obs, clids, object) {
    dists   <- sapply(clids, compute.dist, obs = obs, act.clustering = object)
    return(dists)
    }
  distances <- t(apply(newdata, 1, comp.dists, clids = clids, object = object))
  colnames(distances) <- clids

  # assign clusters with minimal distance   
  assign.cluster    <- function(distances.i) {
    cl <- which.min(distances.i)
    if (length(cl) > 1) cl <- sample(cl, 1) # random sampling in case of non-uniqueness 
    return(cl)
    }
  cls <- apply(distances, 1, assign.cluster)
    
  # create results object
  res <- list()
  res$distances <- distances
  res$cluster   <- cls
  return(res)
  }


########################
### print function
########################

print.orclus <- function(x, ...){
  cat("\nresult from arbitrary ORiented CLUStering\n\n")
  cat("Call:\n")
  print(x$orclus.call)

  cat("\nDimension of subspaces:", x$subspace.dimension, "\n")
  
  cat("\nCluster assignment result: \n")
  print(x$cluster)
  cat("\n")

  cat("\nCorresponding subspaces of the clusters: \n")
  print(x$subspaces)
  cat("\n")

  cat("\nCluster sparsity coefficient: \n")
  print(x$sparsity.coefficient)
  cat("\n\n")
                                                                                    
  invisible(x)
  }

