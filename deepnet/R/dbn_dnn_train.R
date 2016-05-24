##' Training a Deep neural network with weights initialized by DBN
##'
##' Training a Deep neural network with weights initialized by DBN
##' @param x matrix of x values for examples
##' @param y vector or matrix of target values for examples
##' @param hidden vector for number of units of hidden layers.Default is c(10).
##' @param activationfun activation function of hidden unit.Can be "sigm","linear" or "tanh".Default is "sigm" for logistic function
##' @param learningrate learning rate for gradient descent. Default is 0.8.
##' @param momentum momentum for gradient descent. Default is 0.5 .
##' @param learningrate_scale  learning rate will be mutiplied by this scale after every iteration. Default is 1 .
##' @param numepochs number of iteration for samples  Default is 3.
##' @param batchsize size of mini-batch. Default is 100.
##' @param output function of output unit, can be "sigm","linear" or "softmax". Default is "sigm".
##' @param hidden_dropout drop out fraction for hidden layer. Default is 0.
##' @param visible_dropout drop out fraction for input layer Default is 0.
##' @param cd number of iteration for Gibbs sample of CD algorithm.
##' @author Xiao Rong
##' @examples
##' Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' x <- matrix(c(Var1,Var2),nrow=100,ncol=2)
##' y <- c(rep(1,50),rep(0,50))
##' dnn <-dbn.dnn.train(x,y,hidden=c(5,5))
##' ## predict by dnn
##' test_Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' test_Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' test_x <- matrix(c(test_Var1,test_Var2),nrow=100,ncol=2)
##' nn.test(dnn,test_x,y)
##' @export
dbn.dnn.train <- function(x,y,hidden=c(1),
                          activationfun="sigm",
                          learningrate=0.8,
                          momentum=0.5,
                          learningrate_scale=1,
                          output="sigm",
                          numepochs=3,batchsize=100,
                          hidden_dropout=0,visible_dropout=0,cd=1){
  
  output_dim <- 0
  if(is.vector(y)){
    output_dim <- 1
  }else if(is.matrix(y)){
    output_dim <- ncol(y)
  }
  if (output_dim == 0) 
    stop("y must be a vector or matrix!")
  message("begin to train dbn ......")
  dbn <- dbn.train(x,hidden=hidden,
                   numepochs=numepochs,batchsize=batchsize,
                   learningrate=learningrate,learningrate_scale=learningrate_scale,
                   momentum=momentum,cd=cd)
  message("dbn has been trained.")
  initW <- list()
  initB <- list()
  for(i in 1:(length(dbn$size) - 1)){
    initW[[i]] <- dbn$rbm[[i]]$W   
    initB[[i]] <- dbn$rbm[[i]]$C   
  }
  #random init weight between last hidden layer and output layer
  last_hidden <- dbn$size[length(dbn$size)]
  initW[[length(dbn$size)]] <- matrix(runif(output_dim*last_hidden,min=-0.1,max=0.1), c(output_dim,last_hidden))
  initB[[length(dbn$size)]] <- runif(output_dim,min=-0.1,max=0.1)
  message("begin to train deep nn ......")
  dnn <- nn.train(x,y,initW=initW,initB=initB,hidden=hidden,
                  activationfun=activationfun,
                  learningrate=learningrate,
                  momentum=momentum,
                  learningrate_scale=learningrate_scale,
                  output=output,
                  numepochs=numepochs,batchsize=batchsize,
                  hidden_dropout=0,visible_dropout=0)
  message("deep nn has been trained.")
  dnn
}