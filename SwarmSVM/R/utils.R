#' Euclidean Distance calculation
#' 
#' @param x the data matrix
#' @param centers the matrix of centers
#' 
eucliDist = function(x, centers) {
  assertInt(nrow(x), lower = 1)
  assertInt(ncol(x), lower = 1)
  assertMatrix(centers, ncols = ncol(x))
  if (nrow(centers)>1) {
    result = apply(centers, 1, function(C) colSums( (t(x)-C)^2 ))
  } else {
    result = colSums((t(x)-as.vector(centers))^2)
    result = matrix(result, length(result), 1)
  }
  assertMatrix(result, nrows = nrow(x), ncols = nrow(centers))
  return(result)
}

#' Euclidean Distance based clustering prediction
#' 
#' @param x the data matrix
#' @param cluster.object the matrix of centers
#' 
kmeans.predict = function(x, cluster.object) {
  #assertMatrix(x)
  assertInt(nrow(x), lower = 1)
  assertInt(ncol(x), lower = 1)
  centers = cluster.object$centers
  assertMatrix(centers)
  
  euclidean.dist = eucliDist(x, centers)
  result = max.col(-euclidean.dist)
  assertInteger(result, lower = 1, upper = nrow(centers), len = nrow(x))
  return(result)
}

#' Kmeans Clustering from RcppMLPACK
#' 
#' The Kmeans algorithm from \code{RcppMLPACK}.
#' 
#' @param x The input data for the clustering algorithm.
#' @param centers A number indicating the number of clustering centers.
#' @param ... arguments for future use.
#' 
cluster.fun.mlpack.old = function(x, centers, ...) {

  assertInt(nrow(x), lower = 1)
  assertInt(ncol(x), lower = 1)
  assertInt(centers, lower = 1, upper = nrow(x))
  
  BBmisc::suppressAll({
    # result = RcppMLPACK::mlKmeans(t(as.matrix(x)),centers)
  })
  result$cluster = result$result+1
  result$cluster= as.integer(result$cluster)
  k = max(result$cluster)
  cluster.centers = matrix(0,k,ncol(x))
  for (i in 1:k) {
    index = which(result$cluster == i)
    cluster.centers[i,] = colMeans(x[index,,drop = FALSE])
  }
  
  result$centers = cluster.centers
  result$clusters = NULL
  result$result = NULL
  assertMatrix(result$centers, min.rows = 1, ncols = ncol(x))
  assertInteger(result$cluster, lower = 1, upper = nrow(result$centers), 
                len = nrow(x))
  return(result)
}

# Temperal mask of RcppMLPACK
cluster.fun.mlpack = stats::kmeans

cluster.predict.mlpack = function(x, cluster.object) {
  return(kmeans.predict(x, cluster.object))
}

scaleBySD = function(x, sds) {
  assertInt(ncol(x), lower = 1)
  assertInt(nrow(x), lower = 1)
  if (testClass(x, 'dgCMatrix')) {
    x = x %*% Matrix::Diagonal(x = 1/sds)
  } else {
    for (i in 1:ncol(x)) {
      x[,i] = x[,i]/sds[i]
    }
  }
  return(x)
}

sendMsg = function(..., verbose) {
  assertFlag(verbose)
  if (verbose)
    message(...)
}

muteFun = function(expr, mute = FALSE) {
  assertFlag(mute)
  if (mute) {
    BBmisc::suppressAll(expr)
  } else {
    expr
  }
}

#' Wrapper function for kernal kmeans
#' 
#' @param x the input data
#' @param centers the number of centers
#' @param ... other parameters passing to \code{kernlab::kkmeans}
#' 
cluster.fun.kkmeans = function(x, centers, ...) {
  x = as.matrix(x)
  assertMatrix(x, min.rows = 1, min.cols = 1)
  # due to a wierd namespace problem i add this line
  # tmp = kernlab::kkmeans(as.matrix(iris[,-5]), centers, ...)
  assertInt(centers, lower = 1, upper = nrow(x))
  args = list(x = x, centers = centers, ...)
  kernl.result = tryCatch(expr = {kernl.result = do.call(kernlab::kkmeans, args)}, 
                          error = function(err) {
                            missing.flag = grepl('missing', err)
                            compute.flag = grepl('unable', err)
                            err.ind = which(c(missing.flag, compute.flag))
                            if (err.ind == 1) {
                              stop("Too many number of centers, some clusters have zero data point.")
                            } else if (err.ind == 2) {
                              stop("Computational error in kkmeans, please consider changing the seed.")
                            } else {
                              stop(err)
                            }
                          })
  
  result = list()
  result$cluster = kernl.result@.Data
  result$centers = kernl.result@centers
  result$kkmeans.res = kernl.result
  
  return(result)
}

#' Predict function for kernel kmeans
#' 
#' @param x The data to make prediction
#' @param cluster.object The result object from \code{kernlab::kkmeans}
#' 
cluster.predict.kkmeans = function(x, cluster.object) {
  x = as.matrix(x)
  assertMatrix(x, min.rows = 1, min.cols = 1)
  
  kkmeans.res = cluster.object$kkmeans.res
  assertClass(kkmeans.res,'specc')
  
  kern.fun = kkmeans.res@kernelf@.Data
  center.mat = kkmeans.res@centers
  n = nrow(x)
  result = integer(n)
  flag = TRUE
  ind = 1:min(n, 20000)

  while (flag) {
    
    if (max(ind) == n)
      flag = FALSE
    
    kern.mat = kernlab::kernelMatrix(kern.fun, x[ind,,drop = FALSE], center.mat)
    kern.c = diag(kernelMatrix(kern.fun,
                               kkmeans.res@centers,
                               kkmeans.res@centers))
    dist.mat = t(kern.c-2*t(kern.mat))
    result[ind] = max.col(-dist.mat)
    
    ind = (max(ind)+1):min(max(ind)+20000, n)
    
  }
  assertInteger(result, lower = 1, upper = nrow(center.mat), len = nrow(x))
  return(result)
}


#' svmguide1
#' 
#' An astroparticle application from Jan Conrad of Uppsala University, Sweden. 
#' 
#' @docType data
#' @keywords datasets
#' @name svmguide1
#' @usage data(svmguide1)
#' @format A list of two data objects \code{svmguide1} and \code{svmguide1.t}. 
#'    The first column is the target variable.
#' 
NULL


