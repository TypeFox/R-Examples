#' Build a stacked Autoencoder.
#' 
#' @param X.train A matrix of training data.
#' @param n.nodes A vector of numbers containing the number of units 
#'                at each hidden layer.
#' @param unit.type hidden unit activation type as 
#'                per \code{autoencode()} params.
#' @param lambda Vector of scalars indicating weight decay per layer 
#'                as per \code{autoencode()}.
#' @param beta Vector of scalars indicating sparsity penalty per layer
#'                as per \code{autoencode()}.
#' @param rho Vector of scalars indicating sparsity parameter per layer
#'                as per \code{autoencode()}.
#' @param epsilon Vector of scalars indicating initialisation parameter
#'                for weights per layer as per \code{autoencode()}.
#' @param optim.method Optimization method as per \code{optim()}.
#' @param rel.tol Relative convergence tolerance as per \code{optim()}
#' @param max.iterations Maximum iterations for \code{optim()}.
#' @param rescale.flag A logical flag indicating whether input data should
#'                be rescaled.
#' @param rescaling.offset A small non-negative value used for rescaling. 
#'                Further description available in the documentation 
#'                of \code{autoencoder}.
#' @return An object of class \code{SAENET} containing the following elements for each layer of the stacked autoencoder: 
#' \item{ae.out}{ An object of class \code{autoencoder} containing the autoencoder created in that layer of the stacked autoencoder.} 
#' \item{X.output}{ In layers subsequent to the first, a matrix containing the activations of the hidden neurons.}
#' @examples
#' library(autoencoder)
#' data(iris)
#' #### Train a stacked sparse autoencoder with a (5,3) architecture and 
#' #### a relatively minor sparsity penalty. Try experimenting with the 
#' #### lambda and beta parameters if you haven't worked with sparse 
#' #### autoencoders before - it's worth inspecting the final layer
#' #### to ensure that output activations haven't simply converged to the value of 
#' #### rho that you gave (which is the desired activation level on average).
#' #### If the lambda/beta parameters are set high, this is likely to happen.
#' 
#' 
#' output <- SAENET.train(as.matrix(iris[,1:4]), n.nodes = c(5,3), 
#'                        lambda = 1e-5, beta = 1e-5, rho = 0.01, epsilon = 0.01)
#' @export
SAENET.train <- function (X.train, n.nodes = c(4,3,2), unit.type = c("logistic", "tanh"), lambda, beta, rho, epsilon, optim.method = c("BFGS", "L-BFGS-B", "CG"), rel.tol = sqrt(.Machine$double.eps), max.iterations = 2000, rescale.flag = F, rescaling.offset = 0.001){
  h <- as.list(n.nodes)
  if (length(rho)<length(n.nodes)) { rho <- rep(rho,length(n.nodes))}
  if (length(beta)<length(n.nodes)) { beta <- rep(beta,length(n.nodes))}
  if (length(lambda)<length(n.nodes)) { lambda <- rep(lambda,length(n.nodes))}
  if (length(epsilon)<length(n.nodes)) { epsilon <- rep(epsilon,length(n.nodes))}
  for (l in 1:length(n.nodes)){
    if (l == 1){
      ae.out <- autoencoder::autoencode(X.train,N.hidden=n.nodes[l],lambda=lambda[l],epsilon=epsilon[l],beta=beta[l],rho=rho[l],unit.type=unit.type,optim.method = optim.method, rel.tol = rel.tol, max.iterations = max.iterations, rescale.flag = rescale.flag, rescaling.offset = rescaling.offset)
      h[[l]] <- list(ae.out=ae.out,X.output=predict(ae.out,X.train,hidden.output=T)$X.output,n.nodes=n.nodes)
    } else {
      ae.out <- autoencoder::autoencode(h[[l-1]]$X.output,N.hidden=n.nodes[l],lambda=lambda[l],epsilon=epsilon[l],beta=beta[l],rho=rho[l],unit.type=unit.type,optim.method = optim.method, rel.tol = rel.tol, max.iterations = max.iterations, rescale.flag = F)
      h[[l]] <- list(ae.out=ae.out,X.output=predict(ae.out,h[[l-1]]$X.output,hidden.output=T)$X.output)
    }
  }
  class(h) <- "SAENET"
  return(h)
}



#' Use a stacked autoencoder to pre-train a feed-forward neural network.
#' 
#' @param h The object returned from \code{SAENET.train()}
#' @param X.train A matrix of training data.
#' @param target A vector of target values. If given as a factor, classification
#'                will be performed, otherwise if an integer or numeric vector is given
#'                \code{neuralnet()} will default to regression
#' @param nn.control A named list with elements to be passed as control parameters to \code{neuralnet()}. Package defaults used if no values entered.
#' @return An object of class \code{nn} which can be used with the \code{neuralnet} package as normal
#' @export
SAENET.nnet <- function(h,X.train,target,nn.control=NULL){
  if (class(h)!="SAENET"){ message('SAENET.nnet requires a valid object of class SAENET'); break}
  if (is.null(nn.control)) { 
    message('Using defaults supplied by Neuralnet package as no controls detected') 
    neuralnet.control = list(threshold = 0.01, stepmax = 1e+05, rep = 1, startweights = NULL, learningrate.limit = NULL, 
                             learningrate.factor = list(minus = 0.5, plus = 1.2), learningrate = NULL, 
                             lifesign = "none", lifesign.step = 1000, algorithm = "rprop+", 
                             err.fct = "sse", act.fct = "logistic", linear.output = TRUE, 
                             exclude = NULL, constant.weights = NULL, likelihood = FALSE)
  } else {
    neuralnet.control = list(threshold = 0.01, stepmax = 1e+05, rep = 1, startweights = NULL, learningrate.limit = NULL, 
                             learningrate.factor = list(minus = 0.5, plus = 1.2), learningrate = NULL, 
                             lifesign = "none", lifesign.step = 1000, algorithm = "rprop+", 
                             err.fct = "sse", act.fct = "logistic", linear.output = TRUE, 
                             exclude = NULL, constant.weights = NULL, likelihood = FALSE)
    neuralnet.control[names(nn.control)] = nn.control
  }
  #### Prepare weight matrices for neuralnet()
  create.weights <- function(x) {
    nnet.weights <- as.list(rep(0,length(x)))
    for (l in 1:length(x)){
      nnet.weights[[l]] <- t(cbind(x[[l]]$ae.out$b[[1]],x[[l]]$ae.out$W[[1]]))
    }
    return(nnet.weights)
  }
  
  #### Detect what sort of net we're making and reformulate the target variable if necessary
  format.target <- function(x) {
    if (is.factor(x)){
      classes <- levels(x)
      train.target <- as.data.frame(matrix(0,length(x),length(classes)))
      colnames(train.target) <- classes
      for (class in classes){
        train.target[x == class,colnames(train.target) == class] <- 1
      }
      colnames(train.target) <- paste0("target.",classes)
    } else if (is.integer(x) || is.numeric(x)) {
      train.target <- x
    } else {
      message('SAENET detected that the target values are neither integer, factor or numeric')
      break
    }
    return(train.target)
  }
  
  
  nnet.weights <- create.weights(h)
  train.target <- format.target(target)
  output.names <- if (is.data.frame(train.target)){colnames(train.target)}else{"target"}
  #### Pretrain the connections between the last layer of the SAE and the output variables
  X.temp <- h[[length(h)]]$X.output
  X.temp <- as.data.frame(X.temp)
  nnet.weight.train <- cbind(X.temp,train.target)
  train.net <- neuralnet::neuralnet(as.formula(paste(paste(colnames(train.target),collapse=" + "),"~",paste(colnames(X.temp),collapse=" + "))),data=nnet.weight.train,threshold=neuralnet.control$threshold,hidden=0,
                         stepmax=neuralnet.control$stepmax,rep=neuralnet.control$rep,learningrate.limit=neuralnet.control$learningrate.limit,learningrate.factor=neuralnet.control$learningrate.factor,
                         learningrate=neuralnet.control$learningrate,lifesign=neuralnet.control$lifesign,lifesign.step=neuralnet.control$lifesign.step,
                         algorithm=neuralnet.control$algorithm,err.fct=neuralnet.control$err.fct,act.fct=neuralnet.control$act.fct,
                         linear.output=neuralnet.control$linear.output,exclude=neuralnet.control$exclude,constant.weights=neuralnet.control$constant.weights,likelihood=neuralnet.control$likelihood)
  nnet.weights[[length(nnet.weights)+1]] <- train.net$weights[[1]][[1]]
  #### if pre-training the final layer didn't work for some reason, randomly sample instead
  if (length(nnet.weights)==length(h)){
    nnet.weights[[length(nnet.weights)+1]] <- rnorm(length(output.names)*(ncol(nnet.weights[[length(nnet.weights)]])+1),mean=0,sd=0.1^2)
    dim(nnet.weights[[length(nnet.weights)]]) <- c(ncol(nnet.weights[[length(nnet.weights)-1]])+1,length(output.names))
  }
  n.nodes <- unlist(lapply(nnet.weights,function(x) nrow(x)-1)[2:length(nnet.weights)])  
  
  #### Set up the dependent variables
  if (!is.data.frame(X.train)) {X.train <- as.data.frame(X.train)}
  input.data <- cbind(train.target,X.train)
  
  #### Run the neural net
  net <- neuralnet::neuralnet(as.formula(paste(paste(output.names,collapse=" + "),"~",paste(colnames(X.train),collapse=" + "))), data=input.data, startweights=nnet.weights, hidden=n.nodes, threshold=neuralnet.control$threshold,
                   stepmax=neuralnet.control$stepmax,rep=neuralnet.control$rep,learningrate.limit=neuralnet.control$learningrate.limit,learningrate.factor=neuralnet.control$learningrate.factor,
                   learningrate=neuralnet.control$learningrate,lifesign=neuralnet.control$lifesign,lifesign.step=neuralnet.control$lifesign.step,
                   algorithm=neuralnet.control$algorithm,err.fct=neuralnet.control$err.fct,act.fct=neuralnet.control$act.fct,
                   linear.output=neuralnet.control$linear.output,exclude=neuralnet.control$exclude,constant.weights=neuralnet.control$constant.weights,likelihood=neuralnet.control$likelihood)
  return(net)
}


#' Obtain the compressed representation of new data for specified layers from a stacked autoencoder.
#' 
#' @param h The object returned from \code{SAENET.train()}
#' @param new.data A matrix of training data.
#' @param layers A numeric vector indicating which layers of the stacked autoencoder to return output for
#' @param all.layers A boolean value indicating whether to override \code{layers} and return the encoded output for all layers. Defaults to \code{FALSE}
#' @return A list, for which each element corresponds to the output of \code{predict.autoencoder()} from package \code{autoencoder} for the specified layers of the stacked autoencoder. 
#' @examples
#' library(autoencoder)
#' data(iris)
#' #### Train a stacked sparse autoencoder with a (5,3) architecture and 
#' #### a relatively minor sparsity penalty. Try experimenting with the 
#' #### lambda and beta parameters if you haven't worked with sparse 
#' #### autoencoders before - it's worth inspecting the final layer
#' #### to ensure that output activations haven't simply converged to the value of 
#' #### rho that you gave (which is the desired activation level on average).
#' #### If the lambda/beta parameters are set high, this is likely to happen.
#' 
#'
#' output <- SAENET.train(as.matrix(iris[1:100,1:4]), n.nodes = c(5,3), 
#'                        lambda = 1e-5, beta = 1e-5, rho = 0.01, epsilon = 0.01)
#' 
#' 
#' predict.out <- SAENET.predict(output, as.matrix(iris[101:150,1:4]), layers = c(2))
#' @export
SAENET.predict <- function(h, new.data, layers=c(1), all.layers=FALSE){
  if (class(h)!="SAENET"){ message('SAENET.predict requires a valid object of class SAENET'); break}
  if (all.layers==TRUE) {layers <- 1:length(h)}
  results.out <- list(layers)
  for (layer in 1:length(h)){
    if (layer == 1) {
      results.out[[layer]] <- autoencoder::predict.autoencoder(h[[layer]]$ae.out,new.data,hidden.output=TRUE)
    } else if (layer > 1) {
      results.out[[layer]] <- autoencoder::predict.autoencoder(h[[layer]]$ae.out,results.out[[layer - 1]]$X.output,hidden.output=TRUE)
    }
  }
  results.out <- results.out[layers]
  return(results.out)
}
