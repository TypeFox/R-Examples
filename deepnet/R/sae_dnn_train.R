##' Training a Deep neural network with weights initialized by Stacked AutoEncoder
##'
##' Training a Deep neural network with weights initialized by Stacked AutoEncoder
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
##' @param sae_output function of autoencoder output unit, can be "sigm","linear" or "softmax". Default is "linear".
##' @param hidden_dropout drop out fraction for hidden layer. Default is 0.
##' @param visible_dropout drop out fraction for input layer Default is 0.
##' @examples
##' Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' x <- matrix(c(Var1,Var2),nrow=100,ncol=2)
##' y <- c(rep(1,50),rep(0,50))
##' dnn <-sae.dnn.train(x,y,hidden=c(5,5))
##' ## predict by dnn
##' test_Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' test_Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' test_x <- matrix(c(test_Var1,test_Var2),nrow=100,ncol=2)
##' nn.test(dnn,test_x,y)
##' @author Xiao Rong
##' @export
sae.dnn.train <- function(x,y,hidden=c(1),
                          activationfun="sigm",
                          learningrate=0.8,
                          momentum=0.5,
                          learningrate_scale=1,
                          output="sigm",
                          sae_output="linear",
                          numepochs=3,batchsize=100,
                          hidden_dropout=0,visible_dropout=0){
  output_dim <- 0
  if(is.vector(y)){
    output_dim <- 1
  }else if(is.matrix(y)){
    output_dim <- ncol(y)
  }
  if (output_dim == 0) 
    stop("y must be a vector or matrix!")
  message("begin to train sae ......")
  sae <- sae.train(x,hidden=hidden,
                   activationfun=activationfun,
                   output=sae_output,
                   numepochs=numepochs,batchsize=batchsize,
                   learningrate=learningrate,learningrate_scale=learningrate_scale,
                   momentum=momentum,
                   hidden_dropout=hidden_dropout,visible_dropout=visible_dropout)
  message("sae has been trained.")
  initW <- list()
  initB <- list()
  for(i in 1:(length(sae$size) - 1)){
    initW[[i]] <- sae$encoder[[i]]$W[[1]]   
    initB[[i]] <- sae$encoder[[i]]$B[[1]]  
  }
  #random init weights between last hidden layer and output layer
  last_hidden <- sae$size[length(sae$size)]
  initW[[length(sae$size)]] <- matrix(runif(output_dim*last_hidden,min=-0.1,max=0.1), c(output_dim,last_hidden))
  initB[[length(sae$size)]] <- runif(output_dim,min=-0.1,max=0.1)
  message("begin to train deep nn ......")
  dnn <- nn.train(x,y,initW=initW,initB=initB,hidden=hidden,
                  activationfun=activationfun,
                  learningrate=learningrate,
                  momentum=momentum,
                  learningrate_scale=learningrate_scale,
                  output=output,
                  numepochs=numepochs,batchsize=batchsize,
                  hidden_dropout=hidden_dropout,visible_dropout=visible_dropout)
  message("deep nn has been trained.")
  dnn
}