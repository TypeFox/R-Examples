RBFun = function(p,iw,bias){
  #RBF network using Gaussian kernel
  ind = rep(1, dim(p)[1])
  v = matrix(NA, nrow = dim(p)[1], ncol = dim(iw)[1])
  for (i in 1:dim(iw)[1]){
    weight = iw[i,]
    weightMatrix = matrix(rep(weight, length(ind)), ncol = length((weight)), byrow = TRUE)
    v[,i] = -rowSums((p - weightMatrix)^2)
  }
  biasMatrix = matrix(rep(bias, length(ind)), ncol = length((bias)), byrow = TRUE)
  v = v*biasMatrix
  h = exp(v)
}

sinActFun = function(p, iw, bias){
  ### Feedforward neural network using sine activation function ###
  v = as.matrix(p)%*%as.matrix(t(iw))
  ind = rep(1, dim(p)[1])
  biasMatrix = matrix(rep(bias, length(ind)), ncol = length((bias)), byrow = TRUE)
  v = v + biasMatrix
  h = sin(v)
}

sigActFun = function(p, iw, bias){ #toCheck
  ### Feedforward neural network using sigmoidal activation function ###
  v = as.matrix(p)%*%as.matrix(t(iw))
  ind = rep(1, dim(p)[1])
  biasMatrix = matrix(rep(bias, length(ind)), ncol = length((bias)), byrow = TRUE)
  v = v + biasMatrix
  h = 1/(1+exp(-v)); #to check
}

hardLimActFun = function(p, iw, bias){ #to check
  ### Feedforward neural network using hardlim activation function ###
  v = as.matrix(p)%*%as.matrix(t(iw))
  ind = rep(1, dim(p)[1])
  biasMatrix = matrix(rep(bias, length(ind)), ncol = length((bias)), byrow = TRUE)
  v = v + biasMatrix
  h = matrix(as.double(ifelse(v >= 0, 1, 0)), nrow = dim(v)[1], ncol=dim(v)[2])
}

myPinv = function(X, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}

"%^%" <- function(S, power)
  with(eigen(S), vectors %*% (values^power * t(vectors)))

range_01 <- function(x){(x-min(x))/(max(x)-min(x))}
range_11 <- function(x){ 2*(x - min(x))/(max(x) - min(x)) - 1   }

splitData <- function(data, trainPerc = 0.75){
  data = as.data.frame(data)
  data = droplevels(data)
  #set.seed(123)
  smp_size <- floor(trainPerc * nrow(data))
  train_ind <- sample(seq_len(nrow(data)), size = smp_size)
  train <- data[train_ind, ]
  test <- data[-train_ind, ]
  l = list("training" = train, "testing" = test)
  return(l)
}

categoricalToBinary = function(categorical){
  return(model.matrix(~categorical-1))
}

#' Pre processing function for the training and test data set. Each numeric variable is standardized between -1 and 1 and each categorical variable is coded with a dummy coding.
#' @param data to be preprocesses
#' @return return the pre processed dataset
#' @export
preProcess = function(data){
  newData = NULL
  for(i in 1:dim(data)[2]){
    if(is.factor(data[,i])){
      newData = cbind(newData, categoricalToBinary(data[,i]))
    }
    else{
      newData = cbind(newData, range_11(data[,i]))
      colnames(newData)[dim(newData)[2]] = colnames(data)[i]
    }
  }
  return(newData)
}

#' Trains an extreme learning machine with random weights
#' @param formula a symbolic description of the model to be fitted.
#' @param data training data frame containing the variables specified in formula.
#' @param Elm_type select if the ELM must perform a "regression" or "classification"
#' @param nHiddenNeurons number of neurons in the hidden layer
#' @param ActivationFunction "rbf" for radial basis function with Gaussian kernels , "sig" for sigmoidal fucntion, "sin" for sine function, "hardlim" for hard limit function
#' @param N0 size of the first block to be processed
#' @param Block size of each chunk to be processed at each step
#' @return returns all the parameters used in the function, the weight matrix, the labels for the classification, the number of classes found, the bias, the beta activation function and the accuracy on the trainingset
#' @export
OSelm_train.formula = function(formula, data, Elm_type, nHiddenNeurons, ActivationFunction, N0, Block)
{
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  model <- OSelm_training(p = x, y = y, Elm_Type = Elm_type, nHiddenNeurons = nHiddenNeurons, ActivationFunction = ActivationFunction, N0 = N0, Block = Block)
  model$call <- match.call()
  model$formula <- formula
  model
}

