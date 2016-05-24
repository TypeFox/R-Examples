#' Mixture SVMs with gater function
#' 
#' Implementation of Collobert, R., Bengio, S., and Bengio, Y. "A parallel mixture of SVMs for very large scale problems. Neural computation".
#' 
#' @param x the nxp training data matrix. Could be a matrix or an object that can be transformed into a matrix object.
#' @param y a response vector for prediction tasks with one value for each of the n rows of \code{x}. 
#'     For classification, the values correspond to class labels and can be a 1xn matrix, 
#'     a simple vector or a factor. For regression, the values correspond to the values to predict, 
#'     and can be a 1xn matrix or a simple vector.
#' @param m the number of experts
#' @param c a positive constant controlling the upper bound of 
#'     the number of samples in each subset.
#' @param max.iter the number of iterations
#' @param hidden the number of neurons on the hidden layer
#' @param learningrate the learningrate for the back propagation
#' @param threshold neural network stops training once all gradient is below the threshold
#' @param stepmax the maximum iteration of the neural network training process
#' @param seed the random seed. Set it to \code{NULL} to randomize the model.
#' @param valid.x the mxp validation data matrix.
#' @param valid.y if provided, it will be used to calculate the validation score with \code{valid.metric}
#' @param valid.metric the metric function for the validation result. By default it is the accuracy for classification.
#'     Customized metric is acceptable.
#' @param verbose a logical value indicating whether to print information of training.
#' @param ... other parameters passing to \code{neuralnet}
#'
#' @return 
#' \itemize{
#'    \item \code{expert} a list of svm experts
#'    \item \code{gater} the trained neural network model
#'    \item \code{valid.pred} the validation prediction
#'    \item \code{valid.score} the validation score
#'    \item \code{valid.metric} the validation metric
#'    \item \code{time} a list object recording the time consumption for each steps.
#' }
#'
#' @examples
#' 
#' data(svmguide1)
#' svmguide1.t = as.matrix(svmguide1[[2]])
#' svmguide1 = as.matrix(svmguide1[[1]])
#' gaterSVM.model = gaterSVM(x = svmguide1[,-1], y = svmguide1[,1], hidden = 10, seed = 0,
#'                           m = 10, max.iter = 1, learningrate = 0.01, threshold = 1, stepmax = 100,
#'                           valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1], verbose = FALSE)
#' table(gaterSVM.model$valid.pred,svmguide1.t[,1])
#' gaterSVM.model$valid.score
#' 
#' @export
#' 
gaterSVM = function(x, y, m, c = 1, max.iter, 
                    hidden = 5, learningrate = 0.01, threshold = 0.01, stepmax = 100,
                    seed = NULL, valid.x = NULL, valid.y = NULL, valid.metric = NULL,
                    verbose = FALSE, ...) {
  
  total.time.point = proc.time()
  
  if (!testNull(seed))
    set.seed(seed)
  
  assertInt(nrow(x), lower = 1)
  assertInt(ncol(x), lower = 1)
  if (testClass(x, "data.frame"))
    x = data.matrix(x)
  x = as.matrix(x)
  assertClass(x, classes = "matrix")
  n = nrow(x)

  assertVector(y, len = n)
  assertInt(m, lower = 1, upper = n)
  assertInt(c, lower = 0)
  assertInt(max.iter, lower = 1)
  assertInt(hidden, lower = 0)
  assertNumber(learningrate, lower = 1e-8)
  y = as.factor(y)
  if (length(levels(y))!=2)
    stop("Only binary classification is supported")
  y = 2*as.integer(y)-3
  y = matrix(y, length(y), 1)
  # all.data = data.frame(y = as.factor(y), x = x)
  
  # Positive indices
  shuf = sample(which(y == 1))
  fold.len = length(shuf) %/% m
  pos.sub.ind = vector(m, mode = 'list')
  for (i in 1:(m-1)) {
    pos.sub.ind[[i]] = shuf[1:fold.len]
    shuf = setdiff(shuf,pos.sub.ind[[i]])
  }
  pos.sub.ind[[m]] = shuf
  
  # Negative indices
  shuf = sample(which(y == -1))
  fold.len = length(shuf) %/% m
  neg.sub.ind = vector(m, mode = 'list')
  for (i in 1:(m-1)) {
    neg.sub.ind[[i]] = shuf[1:fold.len]
    shuf = setdiff(shuf,neg.sub.ind[[i]])
  }
  neg.sub.ind[[m]] = shuf
  
  # Merge indices
  sub.ind = list()
  for (i in 1:m) {
    sub.ind[[i]] = sample(c(pos.sub.ind[[i]],neg.sub.ind[[i]]))
  }
  
  stopCondition = FALSE
  iter = 1
  S = matrix(0, n, m)
  
  
  svm.time = NULL
  gater.time = NULL
  
  while (!stopCondition) {
    expert = vector(m, mode = 'list')
    constant.pred = rep(0,m)
    time.point = proc.time()
    for (i in 1:m) {
#       sub.data = data.frame(y = as.factor(y[sub.ind[[i]]]), x = x[sub.ind[[i]],])
#       expert[[i]] = alphasvm(y~., data = sub.data, probability = TRUE)
#       S[,i] = predict(expert[[i]], all.data, probability = TRUE)
      #expert[[i]] = alphasvm(x = x[sub.ind[[i]],], y = y[sub.ind[[i]]])
      expert[[i]] = e1071::svm(x = x[sub.ind[[i]],], y = y[sub.ind[[i]]], 
                               scale = FALSE, fitted = FALSE)
      if (expert[[i]]$tot.nSV<1) {
        constant.pred[i] = y[sub.ind[[i]]][1]
        S[,i] = constant.pred[i]
      } else {
        S[,i] = predict(expert[[i]], x)
      }
      #S[,i] = predict(expert[[i]], x)
      sendMsg('Finished training for expert ', i, verbose = verbose)
    }
    svm.time[iter] = (proc.time()-time.point)[3]
    sendMsg('Experts Training Time: ',svm.time[iter], verbose = verbose)
    
    # Train weight
    time.point = proc.time()
    gater.model = gater(x = x, y = y, S = S, hidden = hidden, 
                        learningrate = learningrate, threshold = threshold,
                        verbose = verbose, stepmax = stepmax, ...)
    W = predict(gater.model, x)
    sendMsg('Finish gater training.', verbose = verbose)
    gater.time[iter] = (proc.time()-time.point)[3]
    sendMsg('Gater Neural Net Training Time: ', gater.time[iter], verbose = verbose)
    
    # Re-arrange vectors
    sub.assign = rep(0,n)
    sub.num = rep(0, m)
    sub.avail = rep(1,m)
    W = W-min(W)+1e-9
    for (i in 1:n) {
      ind = which.max(sub.avail*W[i,])
      sub.assign[i] = ind
      sub.num[ind] = sub.num[ind]+1
      if (sub.num[ind]+1 > (n/m+c))
        sub.avail[ind] = -Inf
    }
    for (i in 1:m) {
      sub.ind[[i]] = which(sub.assign==i)
    }
    len = sapply(sub.ind, length)
    if (any(len == 0)) {
      stop("some sub groups have zero length, please adjust the number of experts, or parameter c.")
    }
    
    iter = iter+1
    stopCondition = iter>max.iter
  }
  result = list(expert = expert,
                constant.pred = constant.pred,
                gater = gater.model)
  result = structure(result, class = "gaterSVM")
  
  # Validation
  validation.time = 0
  if (!testNull(valid.x)) {
    time.point = proc.time()
    # assertMatrix(valid.x, min.rows = 1, ncols = ncol(x))
    assertInt(nrow(valid.x), lower = 1)
    assertInt(ncol(valid.x), lower = 1)
    valid.x = as.matrix(valid.x)
    assertClass(valid.x, classes = "matrix")
    if (testNull(valid.y)) {
      warning("Target value for validation is not available.")
      result$valid.pred = predict(result, valid.x)
      result$valid.score = NULL
    } else {
      assertVector(valid.y, len = nrow(valid.x))
      if (testNull((valid.metric))) {
        valid.metric = function(pred, truth) {
          score = sum(diag(table(pred, truth)))/length(truth)
          list(score = score,
               name = 'Accuracy')
        }
      }
      assertFunction(valid.metric)
      result$valid.pred = predict(result, valid.x)
      valid.result = valid.metric(result$valid.pred, valid.y)
      result$valid.score = valid.result$score
      result$valid.metric.name = valid.result$name
    }
    validation.time = (proc.time()-time.point)[3]
    sendMsg("Finished validation process in ", validation.time, ' secs.', 
            verbose = verbose)
  }
  
  total.time = (proc.time()-total.time.point)[3]
  time.list = list(svm.time = svm.time,
                   gater.time = gater.time,
                   validation.time = validation.time,
                   total.time = total.time)
  result$time = time.list

  sendMsg("Finished the whole process in ", total.time, ' secs.', 
          verbose = verbose)
  return(result)
}

#' Prediction for Gater SVM
#' 
#' The function applies a model produced by the 
#'  \code{gaterSVM} function to every row of a data matrix and returns the model predictions.
#' 
#' @param object An object of class \code{gaterSVM}
#' @param newdata newdata An n x p matrix containing the new input data. 
#'     Could be a matrix or a sparse matrix object.
#' @param ... parameters for future usage.
#' 
#' @method predict gaterSVM
#' 
#' @export
#' 
predict.gaterSVM = function(object, newdata, ...) {
  assertClass(object, "gaterSVM")
  assertInt(nrow(newdata), lower = 1)
  if (testClass(newdata, "data.frame"))
    newdata = data.matrix(newdata)
  newdata = as.matrix(newdata)
  
  n = nrow(newdata)
  m = length(object$expert)
  assertInt(m)
  
  S = matrix(0, n, m)
  for (i in 1:m) {
    now.expert = object$expert[[i]]
    if (now.expert$tot.nSV<1) {
      S[,i] = object$constant.pred[i]
    } else {
      S[,i] = predict(now.expert, newdata)
    }
  }
  W = predict(object$gater, newdata)
  pred = sign(object$gater$net$act.fct(rowSums(W*S)))
  return(pred)
}
