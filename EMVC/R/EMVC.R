#
# File: EMVC.R
# Maintainer: rob.frost@dartmouth.edu
# Description: Implementation of the Entropy Minimization over Variable Clusters (EMVC) algorithm.
#              Contains logic for data-driven optimization of 
#              annotations via minimization of the entropy of variable group
#              members over discrete partitions generated via 
#              either k-means clustering or horizontal cuts of dendrograms output via
#              agglomerative hierarchical clustering.
# Copyright (C) Dartmouth College
#

#
# Takes an n-by-p data matrix and a c-by-p binary annotation matrix and generates an optimized, i.e., 
# filtered, version of the annotation matrix by minimizing the entropy between each variable group 
# and the categorical random variable representing membership of each
# variable in the clusters output by either k-means clustering or horizontal cuts of a dendrogram generated via
# agglomerative hierarchical clustering. Annotations are never added during optimization, just removed. 
# Optimization is performed at the level of individual annotations (this has implications for gene sets based on 
# hierarchical ontologies such as GO where parent categories inherit the annotations of child categories).
#
# Args:
#   data:                   Input data matrix, observations-by-variables. Must be specified. Cannot contain missing values.
#   annotations:            Binary category annotation matrix, categories-by-variables. Must be specified.
#   bootstrap.iter:         Number of bootstrap iterations. Defaults to 20. If set to 1, will return the results from 
#                           a single optimization run on the input data matrix (i.e., no bootstrapping will be performed).
#   clust.method            Method used to generate variable clusters. Either "kmeans" or "hclust". Defaults to "kmeans".
#   k.range:                Range of k-means k values or dendrogram cut sizes. Must be specified.
#   hclust.method:          Only relevant if clust.method is "hclust". Method argument for hclust. Defaults to "average".
#   hclust.cor.method       Correlation method used to compute the dissimilarity matrix. Will be supplied as the 
#                           Entries in the dissimilarity matrix will take the form (1-correlation)/2.
#   kmeans.nstart:          Only relevant if clust.method is "kmeans". K-means nstart value. Defaults to 5.
#   kmeans.iter.max:        Only relevant if clust.method is "kmeans".Max number of iterations for k-means. Defaults to 20.
#  
# Returns:
#   optimized.annotations:      Optimized version of the annotation matrix. Contains the average proportion of
#                               cluster sizes in which a given annotation was kept during optimization. 
#                               If bootstrapping is enabled, the optimized 
#                               matrix will contain the average proportions over all bootstrap resampled datasets.
#
EMVC <- function(data, annotations, bootstrap.iter=20, 
    k.range=NA,
    clust.method="kmeans",    
    kmeans.nstart=1,
    kmeans.iter.max=10,
    hclust.method="average",
    hclust.cor.method="spearman") {
  
  params = createEMVCParams(bootstrap.iter=bootstrap.iter,
      k.range=k.range,
      clust.method=clust.method,
      kmeans.nstart=kmeans.nstart,
      kmeans.iter.max=kmeans.iter.max,
      hclust.method=hclust.method)
      
  # Check arguments
  checkOptimizeAnnotationArgs(data, annotations, params)
    
  n = nrow(data)
  p = ncol(data)
  g = nrow(annotations)
  
  # Not bootstrapping, so just directly return results
  if (params$bootstrap.iter <= 1) {
    return (internalOptimize(data, annotations, params))
  }							        

  # Optimize multiple bootstrap resampled datasets
  bootstrap.annotations = matrix(0, nrow=g, ncol=p)
  for (i in 1:params$bootstrap.iter) {    
    if (i %% 10 == 0) {
      message("Bootstrap iteration ", i, ": Sampling ", n, " values with replacement. ",
              "Optimizing ", sum(annotations), " true annotations out of ", (g*p))      
    }
    bootstrap.sample = data[sample(1:n, n, replace=T),]
    opt.annotations = internalOptimize(bootstrap.sample, annotations, params)
    if (i %% 10 == 0) {    
      message("Finished optimization: ", sum(opt.annotations), " annotations out of ", (p*g))
    }
    bootstrap.annotations = bootstrap.annotations + opt.annotations    
  }
    
  # Average the bootstrap proportions
  bootstrap.proportions = bootstrap.annotations/params$bootstrap.iter  
  
  return (bootstrap.proportions)
}

#
# Used to filter the optimized annotation matrix (contains proportions)
# to only contain annotations with a proportion >= the specified value.
#
# Args:
#   annotations:    Optimized annotation matrix with proportions.
#   threshold:      Threshold for annotation filtering.
#
# Returns:
#   Updated version of the input annotations matrix that contains only those annotations
#   with proportions >= the specified threshold.
#
filterAnnotations <- function(annotations, proportion) {
  opt.annotations = apply(annotations, c(1,2), function(x) {
        if (x >= proportion) {
          return (1)
        } 
        return (0)
      })
  return (opt.annotations)	
}

#----------------------------------------------------------------------------------------------------------
# All functions below this line are internal
#----------------------------------------------------------------------------------------------------------

#
# Builds a single list to hold the parameters for the EMVC() method.
#
#
createEMVCParams <- function(bootstrap.iter=20, 
    k.range=NA,
    clust.method="kmeans",    
    kmeans.nstart=1,
    kmeans.iter.max=10,
    hclust.method="average",
    hclust.cor.method="spearman") {
  params = list()
  params$bootstrap.iter = bootstrap.iter
  params$k.range = k.range
  params$clust.method = clust.method
  params$hclust.method = hclust.method
  params$hclust.cor.method = hclust.cor.method  
  params$kmeans.nstart = kmeans.nstart
  params$kmeans.iter.max = kmeans.iter.max  
  return (params)
}

#
# Check arguments for EMVC()
#
checkOptimizeAnnotationArgs <- function(data, annotations, params) {
  current.warn = getOption("warn")
  options(warn=-1)
  if (is.na(data)) {
    stop("data matrix must be specified!")
  }
  if (is.null(annotations)) {
    stop("annotation matrix must be specified!")
  }  
  if (nrow(data) < 2) {
    stop("data matrix must contain at least 2 observations")
  }    
  if (ncol(data) < 4) {
    stop("data matrix must contain at least 4 variables")
  }
  if (ncol(data) < 4) {
    stop("data matrix must contain at least 4 variables")
  }  
  if (ncol(data) != ncol(annotations)) {
    stop("data matrix and annotation matrix must have the same number of columns")
  }
  if (length(which(is.na(data))) > 0) {
    stop("data matrix cannot contain missing values: ", which(is.na(data)))
  }
  if (length(which(is.na(annotations))) > 0) {
    stop("annotation matrix cannot contain missing values")
  }
  if (is.na(params$k.range)) {
    stop("k.range must be specified!")
  }
  if (max(params$k.range) > ncol(data)) {
    stop("max of k.range cannot be larger than the number of variables!")
  }  
  if (min(params$k.range) < 2) {
    stop("min of k.range cannot be less than 2!")
  }  
  if (is.na(params$clust.method)) {
    stop("clust.method must be specified!")
  }
  if (params$clust.method == "kmeans") {
    # all params are optional
  } else if (params$clust.method == "hclust") {
    # all params are optional
  } else {
    stop("clust.method must be set to either 'kmeans' or 'hclust'!")    
  }
  if (params$hclust.cor.method == "spearman") {
  } else if (params$hclust.cor.method == "kendall") {
  } else if (params$hclust.cor.method == "pearson") {    
  } else {
    stop("hclust.cor.method must be set to 'pearson', 'kendall' or 'spearman'!")    
  }  
  
  options(warn=current.warn)  
}

#
# Clusters the variables for the specified data for each k in k.range for the specified cluster method
#
cluster <- function(data, params) {

  #message("Starting clustering...")
  
  p = ncol(data)
  num.cluster.sizes = length(params$k.range)
  
  # will be a matrix with one column per k-value and one row per variable
  # cells hold the cluster number for each variable
  clusters = NA
  
  if (params$clust.method == "kmeans") {
    
    # Execute k-means for each k in k.range
    clusters = matrix(0, nrow=p, ncol=num.cluster.sizes)
    for (i in 1:num.cluster.sizes) {
      k = params$k.range[i]
      #message("Computing k-means for ", k)      
      clusters[,i] = kmeans(t(data), centers=k, nstart=params$kmeans.nstart, iter.max=params$kmeans.iter.max)$cluster
    }
      
  } else if (params$clust.method == "hclust") {

    # Create a correlation-based dissimilarity matrix, using Spearman rank correlation
    cor.diss = as.dist((1-cor(data, method=params$hclust.cor.method))/2)    
    # Cluster variables using hierarchical clustering, correlation dissimilarity and specified method
    hclust.results = hclust(d=cor.diss, method=params$hclust.method)
    # Get all of the partitional cuts
    clusters = cutree(hclust.results, k=params$k.range)  
  }

  #message("...finished clustering")      
  
  return (clusters)
} 


#
# Internal annotation optimization method called during bootstrap optimization
#
internalOptimize <- function(data, annotations, params) {

  n = nrow(data)
  p = ncol(data)
  g = nrow(annotations)  
    
  # Generate partitional clusterings using either k-means or hierarchical clusterings
  clusters = cluster(data, params)
  
  # Create a matrix to hold the sum of optimized annotations
  total.opt.annotations = matrix(0, nrow=g, ncol=p)
  
  # For each cut, compute the entropy-based optimized annotations
  for (i in 1:length(params$k.range)) {
    k = params$k.range[i]    

    # Create the cluster membership matrix, one row per variable, one column per cluster
    cluster.membership = matrix(0, nrow=p, ncol=k)
    for (var in 1:p) {
      if (length(params$k.range) == 1) {
        cluster.membership[var, clusters[var]] = 1
      } else {
        cluster.membership[var, clusters[var,i]] = 1        
      }
    }
    
    # Compute the unoptimized cluster counts. 
    cluster.counts = annotations %*% cluster.membership
        
    # Optimize the annotations based on minimum entropy for each category.
    # This can be simplified to two special cases:
    #
    # 1) One cluster has more annotations than all of the others.
    # 2) Several clusters are tied for the most annotations.
    #
    # For 1): minimum entropy is obtained by simply eliminating all annotations not
    # in that cluster. 
    # For 2): minimum entropy is obtained by doing one of the following:
    #   2a) keep a random selection of the tied clusters and remove all other annotations.
    #   2b) keep the annotations in all tied clusters and remove all other annotations.
    #
    # 2a will actually achieve a minimum entropy solution but the choice is random. 
    # Currently, just 2a is used.
    
    for (j in 1:nrow(cluster.counts)) {
      # Get the index of the cluster with the most annotations for this variable group
      # note: the second random sample is included so that ordering of ties will be random
      largest.cluster = order(cluster.counts[j,], sample(1:k), decreasing=T)[1]
      # Get the indices of the variables in that cluster
      largest.cluster.vars = which((annotations[j,] * cluster.membership[,largest.cluster]) == 1)
      #message("Vars in largest cluster for variable group ", j, ": ", largest.cluster.vars)
      # Add those annotations to the total annotation matrix
      total.opt.annotations[j,largest.cluster.vars] = total.opt.annotations[j,largest.cluster.vars] + 1       
    }    
  }
  
  # Divide the matrix holding the sum of optimized annotations by the number of dendrogram cuts to get the proportion
  # of clusterings in which a given annotation was kept
  opt.annotations = total.opt.annotations/length(params$k.range)
  
  return (opt.annotations)
}
