#' Data Transformation function for Clustered Support Vector Machine
#' 
#' Transform a data matrix according to the kmeans clustering result based on 
#'   Gu, Quanquan, and Jiawei Han. "Clustered support vector machines."
#' 
#' @param x The data matrix, could be a \code{matrix} or \code{dgCMatrix} object.
#' @param lambda The parameter from the algorithm
#' @param cluster.label The clustering label starting from 1. Its length must equal to the number of rows in \code{x}
#' @param sparse Logical argument indicating whether the output should be a sparse matrix or not
#' 
csvmTransform = function(x, lambda, cluster.label, sparse = TRUE) {
  
  assertInt(nrow(x), lower = 1)
  assertInt(ncol(x), lower = 1)
  assertInteger(cluster.label, lower = 1, len = nrow(x))
  assertFlag(sparse)
  
  n = nrow(x)
  m = ncol(x)
  k = max(cluster.label)
  
  if (!(lambda>0))
    stop("Lambda should be strictly larger than zero.")
  
  if (k == 1) {
    warning('Only one cluster for the data, no transform is performed.')
    if (sparse) {
      return(as(x,'matrix.csr'))
    } else {
      return(as.matrix(x))
    }
  }
  
  if (sparse){
    row.index = NULL
    col.index = NULL
    val = NULL
    for (i in 1:k) {
      row.index = c(row.index, rep(which(cluster.label == i),times = m))
      col.index = c(col.index, rep((1 + (i-1)*m):(i*m),each = sum(cluster.label == i)))
      val = c(val, as.vector(x[which(cluster.label == i),]))
    }
    tilde.x = spMatrix(n, k*m, i = row.index, j = col.index, x = val)
    tilde.x = cBind(x / sqrt(lambda), tilde.x)
    tilde.x = as(tilde.x,'dgCMatrix')
    # LiblineaR only support sparse matrix of the class "matrix.csr"
    tilde.x = as(tilde.x, 'matrix.csr')
  } else {
    x = as(x,'matrix')
    tilde.x = matrix(0,n,k*m)
    for (i in 1:k) {
      row.index = which(cluster.label == i)
      col.index = (1 + (i-1)*m):(i*m)
      tilde.x[row.index, col.index] = x[row.index,]
    }
    tilde.x = cbind(x / sqrt(lambda), tilde.x)
  }
  if (sparse) {
    assertClass(tilde.x,"matrix.csr")
  } else {
    assertMatrix(tilde.x)
  }
  return(tilde.x)
}

#' Clustered Support Vector Machine
#' 
#' Implementation of Gu, Quanquan, and Jiawei Han. "Clustered support vector machines."
#' 
#' @param x the nxp training data matrix. Could be a matrix or a sparse matrix object.
#' @param y a response vector for prediction tasks with one value for each of the n rows of \code{x}. 
#'     For classification, the values correspond to class labels and can be a 1xn matrix, 
#'     a simple vector or a factor. For regression, the values correspond to the values to predict, 
#'     and can be a 1xn matrix or a simple vector.
#' @param centers an integer indicating the number of centers in clustering.
#' @param cluster.object an object generated from \code{cluster.fun}, and can be passed to \code{cluster.predict}
#' @param lambda the weight for the global l2-norm
#' @param sparse indicating whether the transformation results in a sparse matrix or not
#' @param valid.x the mxp validation data matrix.
#' @param valid.y if provided, it will be used to calculate the validation score with \code{valid.metric}
#' @param valid.metric the metric function for the validation result. By default it is the accuracy for classification
#'     or RMSE for regression. Customized metric is acceptable.
#' @param type the type of the mission for \code{LiblineaR}.
#' @param cost cost of constraints violation (default: 1). 
#'     Rules the trade-off between regularization and correct classification on data. 
#'     It can be seen as the inverse of a regularization constant. 
#'     See details in \code{LiblineaR}. 
#' @param epsilon set tolerance of termination criterion for optimization. 
#'     If NULL, the LIBLINEAR defaults are used, which are:
#' @param bias if bias is \code{TRUE} (default), instances of data becomes [data; 1].
#' @param wi a named vector of weights for the different classes, 
#'     used for asymmetric class sizes. Not all factor levels have to be supplied (default weight: 1). 
#'     All components have to be named according to the corresponding class label. 
#'     Not used in regression mode.
#' @param verbose if set to 0, no information is printed. 
#'     If set to 1 (default), the running time and validation score (if applicable) will be printed.
#'     If set to 2, the running time ,validation score (if applicable) and the \code{LiblineaR} information will be printed.
#' @param seed the random seed. Set it to \code{NULL} to randomize the model.
#' @param cluster.method The clusterign algorithm to use. Possible choices are 
#' \itemize{
#'     \item "kmeans" Algorithm from \code{stats::kmeans}
#'     \item "mlKmeans" Algorithm from \code{RcppMLPACK::mlKmeans}
#'     \item "kernkmeans" Algorithm from \code{kernlab::kkmeans}
#' }
#' If \code{cluster.fun} and \code{cluster.predict} are provided, \code{cluster.method} doesn't work anymore.
#' @param cluster.fun The function to train cluster labels for the data based on given number of centers. 
#'     Customized function is acceptable, as long as the resulting list contains two fields named as \code{cluster} and \code{centers}.
#' @param cluster.predict The function to predict cluster labels for the data based on trained object. 
#'     Customized function is acceptable, as long as the resulting list contains two fields named as \code{cluster} and \code{centers}.
#' @param ... additional parameters passing to \code{cluster.fun}.
#' 
#' @return 
#' \itemize{
#'    \item \code{svm} the svm object from \code{LiblineaR}
#'    \item \code{lambda} the parameter used.
#'    \item \code{sparse} whether the data is sparsely transformed
#'    \item \code{label} the clustering label for training data
#'    \item \code{centers} the clustering centers from teh training dataset
#'    \item \code{cluster.fun} the function used for clustering
#'    \item \code{cluster.object} the object either
#'    \item \code{cluster.predict} the function used for prediction on new data based on the object
#'    \item \code{valid.pred} the validation prediction
#'    \item \code{valid.score} the validation score
#'    \item \code{valid.metric} the validation metric
#'    \item \code{time} a list object recording the time consumption for each steps.
#' }
#' 
#' @examples
#' data(svmguide1)
#' svmguide1.t = svmguide1[[2]]
#' svmguide1 = svmguide1[[1]]
#' 
#' csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
#'                       centers = 8, seed = 512, verbose = 0,
#'                       valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
#' csvm.pred = csvm.obj$valid.pred
#' 
#' # Or predict from the data
#' pred = predict(csvm.obj, svmguide1.t[,-1])
#' 
#' @export
#' 
clusterSVM = function(x, y, centers = NULL, cluster.object = NULL, lambda = 1, sparse = TRUE, 
                      valid.x = NULL, valid.y = NULL, valid.metric = NULL,
                      type = 1, cost = 1, epsilon = NULL, 
                      bias = TRUE, wi = NULL, verbose = 1, seed = NULL,
                      cluster.method = "kmeans",
                      cluster.fun = NULL, cluster.predict = NULL, ...) {
  
  # Parameter check
  assertInt(nrow(x), lower = 1)
  assertInt(ncol(x), lower = 1)
  if (testClass(x, "data.frame"))
    x = data.matrix(x)
  assertVector(y)
  assertNumber(lambda, lower = 0)
  assertNumber(cost, lower = 0)
  if (!is.null(epsilon)) assertNumber(epsilon, lower = 0)
  assertInt(type, lower = 0, upper = 7)
  assertInt(verbose, lower = 0, upper = 2)
  
  if (testNull(centers) && testNull(cluster.object)) 
    stop('Either number of centers or the clustering result is needed')
  
  if (testNull(cluster.fun) && testNull(cluster.predict)) {
    assertCharacter(cluster.method)
    if (cluster.method == 'kmeans') {
      cluster.fun = stats::kmeans
      cluster.predict = kmeans.predict
    } else if (cluster.method == 'mlKmeans') {
      cluster.fun = cluster.fun.mlpack
      cluster.predict = cluster.predict.mlpack
    } else if (cluster.method == 'kernkmeans') {
      cluster.fun = cluster.fun.kkmeans
      cluster.predict = cluster.predict.kkmeans
    } else {
      stop("Unknow cluster.method.")
    }
  }
  assertFunction(cluster.fun)
  assertFunction(cluster.predict)
  
  if (!is.null(seed))
    set.seed(seed)
  
  # Clustering 
  total.time.point = proc.time()
  time.point = proc.time()
  
  if (testNull(cluster.object)) {
    # Training
    assertInt(centers, lower = 1, upper = nrow(x))
    cluster.result = cluster.fun(x, centers, ...)
    cluster.object = cluster.result
  } else {
    # Predicting
    cluster.result = list()
    cluster.result$cluster = cluster.predict(x, cluster.object, ...)
    cluster.result$centers = cluster.object$centers
  }
  assertMatrix(cluster.result$centers, min.rows = 1, ncols = ncol(x))
  assertInteger(cluster.result$cluster, lower = 1, upper = nrow(cluster.result$centers), 
                len = nrow(x))
  
  cluster.label = cluster.result$cluster
  cluster.centers = cluster.result$centers
  k = nrow(cluster.centers)
  
  clustering.time = (proc.time()-time.point)[3]
  sendMsg('Time for Clustering: ',clustering.time, ' secs', verbose = verbose>0)
  time.point = proc.time()
  
  # Transformation
  tilde.x = csvmTransform(x, lambda, cluster.label, sparse = sparse)
  
  transform.time = (proc.time()-time.point)[3]
  sendMsg('Time for Transforming: ',transform.time, ' secs', verbose = verbose>0)
  time.point = proc.time()
  
  # Training
  # tmp = unique(y)
  svm.result = LiblineaR(data = tilde.x, target = y, type = type, cost = cost, 
                         epsilon = epsilon, bias = bias,
                         wi = wi, cross = 0, verbose = (verbose>=2))
  
  liblinear.time = (proc.time()-time.point)[3]
  sendMsg('Time for Liblinear: ', liblinear.time, ' secs', verbose = verbose>0)
  
  cluster.svm.result = list(svm = svm.result, 
                            lambda = lambda,
                            sparse = sparse,
                            label = cluster.label, 
                            centers = cluster.centers,
                            cluster.fun = cluster.fun,
                            cluster.object = cluster.object,
                            cluster.predict = cluster.predict)
  cluster.svm.result = structure(cluster.svm.result, class = 'clusterSVM')
  
  # Validation
  validation.time = 0
  if (!testNull(valid.x)) {
    time.point = proc.time()
    
    if (testNull(valid.y)) {
      warning("Target value for validation is not available.")
      cluster.svm.result$valid.pred = predict(cluster.svm.result, valid.x)$predictions
      cluster.svm.result$valid.score = NULL
    } else {
      if (testNull(valid.metric)) {
        if (type<=7) {
          valid.metric = function(pred, truth) list(score = sum(pred==truth)/length(truth),
                                                    name = 'Accuracy')
        } else {
          # rmse
          valid.metric = function(pred, truth) list(score = sqrt(mean((pred-truth)^2)),
                                                    name = 'RMSE')
        }
      }
      cluster.svm.result$valid.pred = predict(cluster.svm.result, valid.x)$predictions
      valid.result = valid.metric(cluster.svm.result$valid.pred, valid.y)
      cluster.svm.result$valid.score = valid.result$score
      cluster.svm.result$valid.metric.name = valid.result$name
    }
    
    validation.time = (proc.time()-time.point)[3]
    sendMsg('Time for Validation: ', validation.time, ' secs', verbose = verbose>0)
  }
  
  total.time = (proc.time()-total.time.point)[3]
  sendMsg('\nTotal Time: ', total.time, ' secs\n', verbose = verbose>0)
  if (!is.null(cluster.svm.result$valid.score))
    sendMsg(cluster.svm.result$valid.metric.name, ' Score: ', 
            cluster.svm.result$valid.score, verbose = verbose>0)
  
  time.record = list()
  time.record$clustering.time = clustering.time
  time.record$transform.time = transform.time
  time.record$liblinear.time = liblinear.time
  time.record$validation.time = validation.time
  time.record$total.time = total.time
  
  cluster.svm.result$time = time.record
  assertClass(cluster.svm.result,"clusterSVM")
  return(cluster.svm.result)
}

#' Predictions with Clustered Support Vector Machines
#' 
#' The function applies a model (classification or regression) produced by the 
#'  \code{clusterSVM} function to every row of a data matrix and returns the model predictions.
#' 
#' @param object Object of class "clusterSVM", created by \code{clusterSVM}.
#' @param newdata An n x p matrix containing the new input data. Could be a matrix or a sparse matrix object.
#' @param cluster.predict a function predict new labels on newdata.
#' @param ... other parameters passing to \code{predict.LiblineaR}
#' 
#' @method predict clusterSVM
#' 
#' @export
#' 
predict.clusterSVM = function(object, newdata = NULL, cluster.predict = NULL, ...) {
  
  if (testNull(newdata))
    return(fitted(object$svm))
  
  assertClass(object, 'clusterSVM')
  assertInt(nrow(newdata), lower = 1)
  assertInt(ncol(newdata), lower = 1)
  if (testClass(newdata, "data.frame"))
    newdata = data.matrix(newdata)
  assertMatrix(object$centers, min.rows = 1, ncols = ncol(newdata))
  if (testNull(cluster.predict)) {
    cluster.predict = object$cluster.predict
  }
  assertFunction(cluster.predict)
  assertFlag(object$sparse)
  assertNumber(object$lambda, lower = 0)
  
  # Assign label
  new.labels = cluster.predict(newdata, object$cluster.object)
  
  # Transformation
  tilde.newdata = csvmTransform(newdata, object$lambda, new.labels, object$sparse)
  
  # Make prediction
  preds = predict(object$svm, tilde.newdata, ...)
  return(preds)
}
