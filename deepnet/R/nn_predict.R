##' Predict new samples by Trainded NN
##'
##' Predict new samples by Trainded NN
##' @param nn nerual network trained by function nn.train
##' @param x new samples to predict
##' @return return raw output value of neural network.For classification task,return the probability of a class
##' @examples
##' Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' x <- matrix(c(Var1,Var2),nrow=100,ncol=2)
##' y <- c(rep(1,50),rep(0,50))
##' nn <-nn.train(x,y,hidden=c(5))
##' ## predict by nn
##' test_Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' test_Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' test_x <- matrix(c(test_Var1,test_Var2),nrow=100,ncol=2)
##' yy <- nn.predict(nn,test_x)
##' 
##' @author Xiao Rong
##' @export

nn.predict <- function(nn,x){
  m <- nrow(x)
  post <- x
  #hidden layer
  for(i in 2:(length(nn$size) - 1)){
    pre <- t( nn$W[[i-1]] %*% t(post) + nn$B[[i-1]] )
    if(nn$activationfun == "sigm"){
      post <- sigm( pre )
    }else if(nn$activationfun == "tanh"){
      post <- tanh(pre)
    }else{
      stop("unsupport activation function 'nn$activationfun'");
    }	
    post <- post * (1 - nn$hidden_dropout)
  }
  #output layer
  i <- length(nn$size)
  pre <- t( nn$W[[i-1]] %*% t(post) + nn$B[[i-1]] )
  if(nn$output == "sigm"){
    post <- sigm( pre )
  }else if(nn$output == "linear"){
    post <- pre  
  }else if(nn$output == "softmax"){
    post <- exp(pre)
    post <- post / rowSums(post) 
  }	else{
    stop("unsupport output function!");
  }	
  post
}

##' Test new samples by Trainded NN
##'
##' Test new samples by Trainded NN,return error rate for classification
##' @param nn nerual network trained by function nn.train
##' @param x new samples to predict
##' @param y new samples' label
##' @param t threshold for classification. If nn.predict value >= t then label 1,else label 0
##' @return error rate
##' @examples
##' Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' x <- matrix(c(Var1,Var2),nrow=100,ncol=2)
##' y <- c(rep(1,50),rep(0,50))
##' nn <-nn.train(x,y,hidden=c(5))
##' test_Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' test_Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' test_x <- matrix(c(test_Var1,test_Var2),nrow=100,ncol=2)
##' err <- nn.test(nn,test_x,y)
##' 
##' @author Xiao Rong
##' @export
nn.test <- function (nn,x,y,t=0.5){
  y_p <- nn.predict(nn,x)
  m <- nrow(x)
  y_p[y_p>=t] <- 1
  y_p[y_p<t] <- 0
  error_count <- sum(abs( y_p - y)) / 2
  error_count / m
}